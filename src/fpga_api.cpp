#include "fpga_api.h"
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <cstring>


#define min(x, y) (((x) < (y)) ? (x) : (y))

FPGA::FPGA(off_t data_addr, off_t output_addr, int m_size, int v_size)
{
  m_size_ = m_size;
  v_size_ = v_size;
  data_size_ = (m_size_ + 1) * v_size_; // fpga bram data size

  qvec_ = new char[v_size_];
  qmat_ = new char[m_size_*v_size_];

  m1_size_ = v_size * v_size;
  m2_size_ = v_size * v_size;
  data_size_M = (v_size_+v_size_)*v_size_;
  
  qm1_ = new char[v_size_*v_size_];
  qm2_ = new char[v_size_*v_size_];
  
  qout_ = new int[m_size_];
  qout_M = new int[v_size_*v_size_];

  //output_ = new unsigned int[m_size_]; // use output_ as tempolar output
  output_M = new unsigned int[v_size_*v_size_]; // use output_M as tempolar output

  data_ = new float[data_size_];
  data_M = new float[data_size_M];

  fd_ = open("/dev/mem", O_RDWR);

  qdata_ = new int[data_size_M];
  qdata_M = static_cast<int *>(mmap(NULL, data_size_M, PROT_READ | PROT_WRITE, MAP_SHARED, fd_, data_addr));
  
  output_ = static_cast<unsigned int *>(mmap(NULL, sizeof(unsigned int), PROT_READ | PROT_WRITE, MAP_SHARED, fd_, output_addr));
  num_block_call_ = 0;
}

FPGA::~FPGA()
{
  munmap(qdata_M, data_size_);
  munmap(output_, sizeof(unsigned int));
  close(fd_);

  delete[] output_;
  delete[] data_;
  delete[] qvec_;
  delete[] qmat_;
  delete[] qout_;
}

float *FPGA::matrix(void)
{
  return data_ + v_size_;
}

float *FPGA::vector(void)
{
  return data_;
}

float *FPGA::matrix_M1(void)
{
  return data_M;
}

float *FPGA::matrix_M2(void)
{
  return data_M + m1_size_;
}

int *FPGA::qmatrix(void)
{
  return qdata_ + v_size_;
}

int *FPGA::qvector(void)
{
  return qdata_;
}

int *FPGA::qmatrix_M1(void)
{
  return qdata_M;
}

int *FPGA::qmatrix_M2(void)
{
  return qdata_M + m1_size_;
}

void FPGA::reset(void)
{
  num_block_call_ = 0;
}

int FPGA::num_block_call(void)
{
  return num_block_call_;
}

void quantize(const float* input, char* quantized, int num_input, int bits_min, int bits_max, int offset, float scale)
{
  for(int i = 0; i < num_input; i++)
  {
    // ceiling function from cmath library (standard lib).
    int quant = ceil(input[i] / scale) + offset;
    if(quant > bits_max)
      quant = bits_max;
    else if(quant < bits_min)
      quant = bits_min;
    
    quantized[i] = quant - offset; // TODO: convert floating point to quantized value
  }
}

void dequantize(const int* quantized, float* output, int num_output, int offset, float scale)
{
  for(int i = 0; i < num_output; i++)
  {
    output[i] = scale * ((float) (quantized[i] - offset)); // TODO: convert quantized value to floating point
  }
}

const int *__attribute__((optimize("O0"))) FPGA::qblockMM(Compute* comp)
{
  num_block_call_ += 1;

  float* m1 = this->matrix_M1();
  float* m2 = this->matrix_M2();
  char* qm1 = reinterpret_cast<char*>(qm1_);
  char* qm2 = reinterpret_cast<char*>(qm2_);
  int* qm1t = this->qmatrix_M1();
  int* qm2t = this->qmatrix_M2();
  float act_scale = 0;
  float weight_scale = 0;
  memset(qdata_, 0, sizeof(int) * data_size_M);

  if(comp->quantized)
  {
    char act_bits_min = 0;
    char act_bits_max = (1<<(comp->act_bits-1))-1;

    float act_min = m2[0];
    float act_max = m2[0];
    for(int i=0; i<v_size_*v_size_; i++){
      if(m2[i] < act_min)
        act_min = m2[i];
      else if(m2[i] > act_max)
        act_max = m2[i];
    }

    act_scale = (act_max - act_min) / 127; // TODO calculate the scale factor
    char act_offset = (char) ceil(-act_min / act_scale); // TODO calculate the zero-offset
    quantize(m2, qm2, v_size_*v_size_, act_bits_min, act_bits_max, act_offset, act_scale); // TODO complete quantize function

    char weight_bits_min = 0;
    char weight_bits_max = (1<<(comp->weight_bits-1))-1;

    float weight_min = m1[0];
    float weight_max = m1[0];
    for(int i=0; i<v_size_*v_size_; i++){
      if(m1[i] < weight_min)
        weight_min = m1[i];
      else if(m1[i] > weight_max)
        weight_max = m1[i];
    }
    weight_scale = (weight_max - weight_min) / 127; // TODO calculate the scale factor
    char weight_offset = (char) ceil(-weight_min / weight_scale); // TODO calculate the zero-offset
    quantize(m1, qm1, v_size_*v_size_, weight_bits_min, weight_bits_max, weight_offset, weight_scale); // TODO complete quantize function

    for(int x=0; x < v_size_/8; x++)
      for(int y=0; y < v_size_/8; y++)
        for(int z=0; z < v_size_/8; z++){
          for(int i=0; i<8; i++)
            for(int j=0; j<8; j++){
              qm1t[i*8 + j] = (int) qm1[(x*8+i)*v_size_ + z*8 + j];
              qm1t[i*8 + j + 64] = (int) qm2[(z*8+i)*v_size_ + y*8 + j];
            }
         
          *output_ = 0x5555;
          while (*output_ == 0x5555);

          for(int i=0; i<8; i++)
            for(int j=0; j<8; j++){
              qdata_[(x*8+i)*v_size_ + y*8 + j] += qdata_M[i*8+j];
            }
        }
    }

    data_M[0] = act_scale * weight_scale;
    return qdata_;
}

const float* FPGA::blockMM(Compute* comp)
{
  num_block_call_ += 1;

  // cpu version
  float* m1 = this->matrix_M1();
  float* m2 = this->matrix_M2();
  char* qm1 = reinterpret_cast<char*>(qm1_);
  char* qm2 = reinterpret_cast<char*>(qm2_);
  float* out  = reinterpret_cast<float*>(output_M);
  // char* qout_M = reinterpret_cast<char*>(qout_M);

  if(comp->quantized)
  {
    char act_bits_min = 0;
    char act_bits_max = (1<<(comp->act_bits-1))-1;

    float act_min = m2[0];
    float act_max = m2[0];
    for(int i=0; i<v_size_*v_size_; i++){
      if(m2[i] < act_min)
        act_min = m2[i];
      else if(m2[i] > act_max)
        act_max = m2[i];
    }
    float act_scale = (act_max - act_min) / 127; // TODO calculate the scale factor
    char act_offset = (char) ceil(-act_min / act_scale); // TODO calculate the zero-offset
    quantize(m2, qm2, v_size_*v_size_, act_bits_min, act_bits_max, act_offset, act_scale); // TODO complete quantize function

    char weight_bits_min = 0;
    char weight_bits_max = (1<<(comp->weight_bits-1))-1;

    float weight_min = m1[0];
    float weight_max = m1[0];
    for(int i=0; i<v_size_*v_size_; i++){
      if(m1[i] < weight_min)
        weight_min = m1[i];
      else if(m1[i] > weight_max)
        weight_max = m1[i];
    }
    float weight_scale = (weight_max - weight_min) / 127; // TODO calculate the scale factor
    char weight_offset = (char) ceil(-weight_min / weight_scale); // TODO calculate the zero-offset
    quantize(m1, qm1, v_size_*v_size_, weight_bits_min, weight_bits_max, weight_offset, weight_scale); // TODO complete quantize function

    for(int i = 0; i < v_size_; ++i)
    {
      for(int j = 0; j < v_size_; ++j){
        qout_M[v_size_*i+j] = 0;
        for(int k = 0; k < v_size_; ++k){
          qout_M[v_size_*i+j] += (qm1[v_size_*i+k]) * (qm2[v_size_*k + j]);
        }
      }
    }

    float output_scale = act_scale * weight_scale;
    dequantize(qout_M, out, v_size_*v_size_, 0, output_scale); // TODO complete dequantize function

  }
  else{
    for(int i = 0; i < v_size_; ++i)
    {
      for(int j = 0; j < v_size_; ++j){
        out[v_size_*i+j] = 0;
        for(int k = 0; k < v_size_; ++k){
          out[v_size_*i+j] += m1[v_size_*i+k] * m2[v_size_*k + j];
        }
      }
    }
  }

  for(int i = 0; i < m1_size_; ++i)
    data_M[i] = out[i];

  return data_M;
}

const float *FPGA::blockMV(Compute* comp)
{
  num_block_call_ += 1;

  // cpu version
  float *vec = this->vector();
  float *mat = this->matrix();
  float *out = reinterpret_cast<float *>(qout_);

  if(comp->quantized)
  {
    char act_bits_min = 0;
    char act_bits_max = (1<<(comp->act_bits-1))-1;

    float act_min = vec[0];
    float act_max = vec[0];
    for(int i=0; i<v_size_; i++){
      if(vec[i] < act_min)
        act_min = vec[i];
      else if(vec[i] > act_max)
        act_max = vec[i];
    }
    float act_scale = (act_max - act_min) / 127; // TODO calculate the scale factor
    char act_offset = (char) ceil(-act_min / act_scale); // TODO calculate the zero-offset
    quantize(vec, qvec_, v_size_, act_bits_min, act_bits_max, act_offset, act_scale); // TODO complete quantize function

    char weight_bits_min = 0;
    char weight_bits_max = (1<<(comp->weight_bits-1))-1;

    float weight_min = mat[0];
    float weight_max = mat[0];
    for(int i=0; i<m_size_*v_size_; i++){
      if(mat[i] < weight_min)
        weight_min = mat[i];
      else if(mat[i] > weight_max)
        weight_max = mat[i];
    }
    float weight_scale = (weight_max - weight_min) / 127; // TODO calculate the scale factor
    char weight_offset = (char) ceil(-weight_min / weight_scale); // TODO calculate the zero-offset
    quantize(mat, qmat_, v_size_*v_size_, weight_bits_min, weight_bits_max, weight_offset, weight_scale); // TODO complete quantize function

    for (int i = 0; i < m_size_; ++i)
    {
      qout_[i] = 0;
      for (int j = 0; j < v_size_; ++j)
        qout_[i] += (qvec_[j]) * (qmat_[v_size_ * i + j]);
    }

    dequantize(qout_, out, m_size_, 0, act_scale * weight_scale); // TODO complete dequantize function
  }
  else
  {
    for (int i = 0; i < m_size_; ++i)
    {
      out[i] = 0;
      for (int j = 0; j < v_size_; ++j)
        out[i] += vec[j] * mat[v_size_ * i + j];
    }
  }

  for (int i = 0; i < m_size_; ++i)
    data_[i] = out[i];

  return data_;
}

void FPGA::largeMM(const float* weight_mat, const float* input_mat, float* output, int num_input, int num_output, int num_matrix2, Compute* comp)
{
  float* m1 = this->matrix_M1();
  float* m2 = this->matrix_M2();

  // 0) Initialize output vector		
  for(int i = 0; i < num_output*num_matrix2; ++i)
    output[i] = 0;

  for(int i = 0; i < num_output; i += v_size_)
  {
    for(int j = 0; j < num_input; j += v_size_)
    {			
      for(int k = 0; k < num_matrix2; k += v_size_)
      {
        // 0) Initialize input vector
        int block_row = min(v_size_, num_output-i);
        int block_col_1 = min(v_size_, num_input-j);
        int block_col_2 = min(v_size_, num_matrix2-k);

        // 1) Assign a m1  
        memset(m1, 0, sizeof(float) * v_size_ * v_size_);
        for(int x=0; x<block_row; x++){
          for(int y=0; y<block_col_1; y++)
            m1[x*v_size_ + y] = weight_mat[(i+x)*num_input + j + y];
        }

        // 2) Assign a m2 
        memset(m2, 0, sizeof(float) * v_size_ * v_size_);
        for(int x=0; x<block_col_1; x++){
          for(int y=0; y<block_col_2; y++)
            m2[x*v_size_ + y] = input_mat[(j+x)*num_matrix2 + k + y];
        }

        // 3) Call a function `blockMM() to execute Matrix matrix multiplication
        const int* ret = this->qblockMM(comp);

        float output_scale = data_M[0];
        float* rett = new float[v_size_*v_size_];
        dequantize(ret, rett, v_size_*v_size_, 0, output_scale);
      
        // 4) Accumulate intermediate results
        for(int n = 0; n<block_row; ++n)
        {
          for(int m = 0; m<block_col_2; ++m)
          {
            float ret_tmp = rett[n*v_size_ + m];
            if(ret_tmp > 10 || ret_tmp < -10)
              ret_tmp = 0;

            output[(i + n) + (k + m)*num_output] += ret_tmp;
          }
        }
      }
    } 
  }
}

void FPGA::largeMV(const float *large_mat, const float *input, float *output, int num_input, int num_output, Compute* comp)
{
  float *vec = this->vector();
  float *mat = this->matrix();

  // 0) Initialize output vector
  for (int i = 0; i < num_output; ++i)
    output[i] = 0;

  for (int i = 0; i < num_output; i += m_size_)
  {
    for (int j = 0; j < num_input; j += v_size_)
    {
      // 0) Initialize input vector
      int block_row = min(m_size_, num_output - i);
      int block_col = min(v_size_, num_input - j);

      // 1) Assign a vector
      memset(vec, 0, sizeof(float) * v_size_);
      for(int k=0; k<block_col; k++)
        vec[k] = input[j+k];

      // 2) Assign a matrix 
      memset(mat, 0, sizeof(float) * v_size_ * m_size_);
      for(int k=0; k<block_row; k++)
        for(int l=0; l<block_col; l++)
          mat[k*v_size_ + l] = large_mat[(i+k)*num_input + j + l];

      // 3) Call a function `blockMV() to execute MV multiplication
      const float* ret = this->blockMV(comp);

      // 4) Accumulate intermediate results
      for (int row = 0; row < block_row; ++row)
        output[i + row] += ret[row];
    }
  }
}

void FPGA::convLowering(const std::vector<std::vector<std::vector<std::vector<float>>>> &cnn_weights,
                        std::vector<std::vector<float>> &new_weights,
                        const std::vector<std::vector<std::vector<float>>> &inputs,
                        std::vector<std::vector<float>> &new_inputs)
{
  /*
   * Arguments:
   *
   * conv_weights: [conv_channel, input_channel, conv_height, conv_width]
   * new_weights: [?, ?]
   * inputs: [input_channel, input_height, input_width]
   * new_inputs: [?, ?]
   *
   */

  int conv_channel = cnn_weights.size();
  int input_channel = cnn_weights[0].size();
  int conv_height = cnn_weights[0][0].size();
  int conv_width = cnn_weights[0][0][0].size();
  int input_height = inputs[0].size();
  int input_width = inputs[0][0].size();
  int num_h = input_height - conv_height + 1;
  int num_w = input_width - conv_width + 1;

  // IMPLEMENT THIS
  // Assign conv_weights
  for(int c1=0; c1<conv_channel; c1++)
    for(int c2=0; c2<input_channel; c2++)
      for(int i=0; i<conv_height; i++)
        for(int j=0; j<conv_width; j++)
          new_weights[c1][c2*conv_height*conv_width + i*conv_width + j] = cnn_weights[c1][c2][i][j];

  // Assign new_inputs
  for(int c=0; c<input_channel; c++)
    for(int i=0; i<num_h; i++)
      for(int j=0; j<num_w; j++)
        for(int x=0; x<conv_height; x++)
          for(int y=0; y<conv_width; y++)
            new_inputs[c*conv_height*conv_width + x*conv_width + y][i*num_w + j] = inputs[c][i+x][j+y];
}
