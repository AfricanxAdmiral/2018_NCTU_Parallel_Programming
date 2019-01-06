#include <cstdio>
#include <cmath>
#include <complex>
#include <cstring>
#include <iostream>
#include <fstream>
using namespace std;

const int N = (1 << 30);

__global__ void multiply(int n, int m, char x[], char y[], int ans[]) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if( i < n && j < m){
        int a = n-1-i, b = m-1-j;
        atomicAdd(&ans[a+b], (x[i]-48) * (y[j]-48));
        //printf("x[%d] : %d * y[%d] : %d = %d\n", i, x[i]-48, j, y[j]-48, ans[a+b]);
    }
}

int main(int argc,char *argv[]) {
    char *s = new char[N];
    char *t = new char[N];

    int n, m, l;

    FILE* fin = fopen(argv[1],"r");
    fscanf(fin,"%s",s);
    fclose(fin);

    fin = fopen(argv[2],"r");
    fscanf(fin,"%s",t);
    fclose(fin);

    n = strlen(s);
    m = strlen(t);
    l = n + m + 1;

    int *ans = new int[l];

    char* cuda_s;
    char* cuda_t;
    int* cuda_ans;
    cudaMalloc(&cuda_s, n * sizeof(char));
    cudaMalloc(&cuda_t, m * sizeof(char));
    cudaMalloc(&cuda_ans, l * sizeof(int));

    cudaMemcpy(cuda_s, s, n*sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_t, t, m*sizeof(char), cudaMemcpyHostToDevice);

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks(n+16 / threadsPerBlock.x, m+16 / threadsPerBlock.y);
    multiply<<<numBlocks, threadsPerBlock>>>(n, m, cuda_s, cuda_t, cuda_ans);

    cudaMemcpy(ans, cuda_ans, l*sizeof(int), cudaMemcpyDeviceToHost);

    for (int i = 0; i < l; ++i) {
        ans[i + 1] += ans[i] / 10;
        ans[i] %= 10;
    }

    int p = l;

    for (; p && !ans[p]; --p);
    for (; ~p; putchar(ans[p--] + '0'));
    puts("");

    return 0;
}
