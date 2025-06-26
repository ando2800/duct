typedef size_t ptr_t;
// not used 
#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */
  int transpose_gpu(ptr_t *handle, int *Nx, int *Ny, const double *alpha, const ptr_t *devPtrA, ptr_t *devPtrC);
#if defined(__cplusplus)
}
#endif /* __cplusplus */
