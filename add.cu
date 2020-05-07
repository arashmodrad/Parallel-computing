//Just your regular Hello World file
// to be compiled with nvcc rather than gcc

#include <stdio.h>

void add_on_CPU(float *a, float *b, float *c){
  *c = *a + *b;
}

__global__ void add_on_GPU(float *a, float *b, float *c) {
  *c = *a + *b;
}


int main(void) {
  float a, b, c;            // host copies of a, b, c
  float *d_a, *d_b, *d_c;     // device copies of a, b, c
  int size = sizeof(float);
  
  // Allocate space for device copies of a, b, c
  cudaMalloc((float **)&d_a, size);
  cudaMalloc((float **)&d_b, size);
  cudaMalloc((float **)&d_c, size);
  
  // Setup input values
  a = 2.1;
  b = 7;

  //compute the sum on CPU
  add_on_CPU(&a,&b,&c);
  printf("Result on CPU = %g\n",c);

  // Copy inputs to device
  cudaMemcpy(d_a, &a, size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_b, &b, size, cudaMemcpyHostToDevice);
  
  // Launch add_on_GPU() kernel on GPU
  add_on_GPU<<<1,1>>>(d_a, d_b, d_c);
  
  // Copy result back to host
  cudaMemcpy(&c, d_c, size, cudaMemcpyDeviceToHost);

  printf("Result on GPU = %g\n",c);

  // Cleanup
  cudaFree(d_a); cudaFree(d_b); cudaFree(d_c);
  
  return 0;
}
