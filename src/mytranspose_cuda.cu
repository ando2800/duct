#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#if defined(__GNUC__)
#include <stdint.h>
#endif /* __GNUC__ */
#include "cuda.h" 
#include "cuda_runtime.h" 

//#include <cublas.h>
#include <cublas_v2.h>

//#include "transpose_cuda_common.h"
//#define TRANSPOSE_GPU             transpose_gpu_
//#define TRANSPOSE_GPU             TRANSPOSE_GPU
//#include "transpose_cuda.h" // not necessary 

#define TILE_DIM 32
#define BLOCK_ROWS 8

typedef size_t ptr_t;

__global__ void copy(double *odata, const double *idata)
{
  int x = blockIdx.x * TILE_DIM + threadIdx.x;
  int y = blockIdx.y * TILE_DIM + threadIdx.y;
  int width = gridDim.x * TILE_DIM;

  for (int j = 0; j < TILE_DIM; j+= BLOCK_ROWS)
    odata[(y+j)*width + x] = idata[(y+j)*width + x];
}

__global__ void matmul(double *odata, const double *idata, const double *idata2)
{
  int x = blockIdx.x * TILE_DIM + threadIdx.x;
  int y = blockIdx.y * TILE_DIM + threadIdx.y;
  int width = gridDim.x * TILE_DIM;

  for (int j = 0; j < TILE_DIM; j+= BLOCK_ROWS)
    odata[(y+j)*width + x] = idata[(y+j)*width + x] * idata2[(y+j)*width + x];
}

__global__ void mat_invdiag(double *odata, const double *F, const double *ED2, const double xlam)
{
  int x = blockIdx.x * TILE_DIM + threadIdx.x;
  int y = blockIdx.y * TILE_DIM + threadIdx.y;
  int width = gridDim.x * TILE_DIM;

  for (int j = 0; j < TILE_DIM; j+= BLOCK_ROWS)
    odata[(y+j)*width + x] = F[(y+j)*width + x] / ( ED2[(y+j)*width + x] - xlam );
}



__global__ void transposeNaive(double *odata, double * idata,
                               int width, int height)
{
  int xIndex = blockIdx.x*TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y*TILE_DIM + threadIdx.y;
  int index_in = xIndex + width * yIndex;
  int index_out = yIndex + height * xIndex;
  for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
    odata[index_out+i] = idata[index_in+i*width];
  }
}

__global__ void transposeCoalesced(double *odata,
                                   double *idata, int width, int height)
{
  __shared__ double tile[TILE_DIM][TILE_DIM];
  int xIndex = blockIdx.x*TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y*TILE_DIM + threadIdx.y;
  int index_in = xIndex + (yIndex)*width;
  xIndex = blockIdx.y * TILE_DIM + threadIdx.x;
  yIndex = blockIdx.x * TILE_DIM + threadIdx.y;
  int index_out = xIndex + (yIndex)*height;
  for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
    tile[threadIdx.y+i][threadIdx.x] =
      idata[index_in+i*width];
  }
  __syncthreads();
  for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
    odata[index_out+i*height] =
      tile[threadIdx.x][threadIdx.y+i];
  }
}

__global__ void transposeFineGrained(double *odata,
                                     double *idata, int width, int height)
{
  __shared__ double block[TILE_DIM][TILE_DIM+1];
  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  int index = xIndex + (yIndex)*width;
  for (int i=0; i < TILE_DIM; i += BLOCK_ROWS) {
    block[threadIdx.y+i][threadIdx.x] =
      idata[index+i*width];
  }
  __syncthreads();
  for (int i=0; i < TILE_DIM; i += BLOCK_ROWS) {
    odata[index+i*height] =
      block[threadIdx.x][threadIdx.y+i];
  }
}

__global__ void transposeCoarseGrained(double *odata,
                                       double *idata, int width, int height)
{
  __shared__ double block[TILE_DIM][TILE_DIM+1];
  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  int index_in = xIndex + (yIndex)*width;
  xIndex = blockIdx.y * TILE_DIM + threadIdx.x;
  yIndex = blockIdx.x * TILE_DIM + threadIdx.y;
  int index_out = xIndex + (yIndex)*height;
  for (int i=0; i<TILE_DIM; i += BLOCK_ROWS) {
    block[threadIdx.y+i][threadIdx.x] =
      idata[index_in+i*width];
  }
  __syncthreads();
  for (int i=0; i<TILE_DIM; i += BLOCK_ROWS) {
    odata[index_out+i*height] =
      block[threadIdx.y+i][threadIdx.x];
  }
}

// interfaces 
extern "C" {
  //int TRANSPOSE_GPU(ptr_t *handle, const int *Nx, const int *Ny, const double *alpha, const ptr_t *devPtrA, ptr_t *devPtrC)
  int transpose_gpu_(ptr_t *handle, const int *Nx, const int *Ny, const double *alpha, const ptr_t *devPtrA, ptr_t *devPtrC)
  {
    
    double *A = (double *)(*devPtrA);
    double *C = (double *)(*devPtrC);
    const int size_x = *Nx;
    const int size_y = *Ny;
    
    // execution configuration parameters
    dim3 grid(size_x/TILE_DIM, size_y/TILE_DIM),
      threads(TILE_DIM,BLOCK_ROWS);
    
    // copy<<<grid, threads>>>(C, A); // OK
    // transposeNaive<<<grid, threads>>>(C, A, size_x, size_y); // OK on vermeer
    transposeCoalesced<<<grid, threads>>>(C, A, size_x, size_y); // OK on vermeer
    // transposeCoarseGrained<<<grid, threads>>>(C, A, size_x, size_y); // testing the performance, result is not correct 
    // transposeFineGrained<<<grid, threads>>>(C, A, size_x, size_y); // testing the performance, result is not correct 
    return 0;
    
    // I could not call cubasDgeam successfully, which returns error code (1) always.
    //  return (int)cublasDgeam((cublasHandle_t)(*handle),CUBLAS_OP_T,CUBLAS_OP_T,*Nx,*Ny,
    //                    alpha,A, *Nx,
    //                    0, 0, *Nx,
    //                    C,*Nx); 
    
    //if (stat != CUBLAS_STATUS_SUCCESS) {
    //  printf ("CUBLAS DGEAM failed\n");
    //return EXIT_FAILURE;
    //} 
    
  }

  
  int matmul_gpu_(ptr_t *handle, const int *Nx, const int *Ny, const double *alpha, const ptr_t *devPtrA, const ptr_t *devPtrB, ptr_t *devPtrC)
  {
    
    double *A = (double *)(*devPtrA);
    double *B = (double *)(*devPtrB);
    double *C = (double *)(*devPtrC);
    const int size_x = *Nx;
    const int size_y = *Ny;
    
    // execution configuration parameters
    dim3 grid(size_x/TILE_DIM, size_y/TILE_DIM),
      threads(TILE_DIM,BLOCK_ROWS);
    matmul<<<grid, threads>>>(C, A, B); // OK on vermeer
    return 0;
  }

  int mat_invdiag_gpu_(ptr_t *handle, const int *Nx, const int *Ny, const double *xlam, const ptr_t *devPtrF, const ptr_t *devPtrED, ptr_t *devPtrC)
  {
    
    double *F = (double *)(*devPtrF);
    double *ED = (double *)(*devPtrED);
    double *C = (double *)(*devPtrC);
    const int size_x = *Nx;
    const int size_y = *Ny;
    const double lam = *xlam;

    // execution configuration parameters
    dim3 grid(size_x/TILE_DIM, size_y/TILE_DIM),
      threads(TILE_DIM,BLOCK_ROWS);
    mat_invdiag<<<grid, threads>>>(C, F, ED, lam); 
    return 0;
  }

  int cuda_stream_create_(cudaStream_t **stream)
  {
    *stream = (cudaStream_t *) malloc(sizeof(cudaStream_t));
    return cudaStreamCreate(*stream);
  }
  
  int cublas_set_stream_(cublasHandle_t *handle, cudaStream_t *streamid)
  {
    return cublasSetStream(*handle, *streamid);
  }
  
  void cuda_stream_destroy_(cudaStream_t *stream)
  {
    cudaStreamDestroy(*stream);
  }
  
}
