#include <cstdio>
#include <cmath>
#include <complex>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <fstream>
using namespace std;

void multiply(int n, int m, int x[], int y[], int ans[]) {
    for (int i = n-1; i >= 0; i--) {
        for (int j = m-1;j >= 0; j--) {
            int a = n-1-i, b = m-1-j;
            ans[a+b] = ans[a+b] + x[i] * y[j];            
        }
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
    l = trans(n + m - 1);

    char *ans = new char[l];

    multiply(n, m, s, t, ans);

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