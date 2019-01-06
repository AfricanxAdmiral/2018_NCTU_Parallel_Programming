#include <cstdio>
#include <cmath>
#include <complex>
#include <cstring>
#include<iostream>
#include <omp.h>
#include<fstream>
using namespace std;
const double PI(acos(-1.0));
typedef complex<double> C;

const int N = (1 << 30);



int trans(int x) {
  return 1 << int(ceil(log(x) / log(2) - 1e-9));  // math.h/log() 以e为底
}
 double sgn;
  double *w;
  double wtime;
  double *x,*x1;
  double *y,*y2;
  double *z;
  double z0;
  double z1;
void cffti ( int n, double w[] )
{
  double arg;
  double aw;
  int i;
  int n2;
  const double pi = 3.141592653589793;

  n2 = n / 2;
  aw = 2.0 * pi / ( ( double ) n );

# pragma omp parallel \
    shared ( aw, n, w ) \
    private ( arg, i )

# pragma omp for nowait

  for ( i = 0; i < n2; i++ )
  {
    arg = aw * ( ( double ) i );
    w[i*2+0] = cos ( arg );
    w[i*2+1] = sin ( arg );
  }
  return;
}


void step ( int n, int mj, double a[], double b[], double c[],
  double d[], double w[], double sgn )
{
  double ambr;
  double ambu;
  int j;
  int ja;
  int jb;
  int jc;
  int jd;
  int jw;
  int k;
  int lj;
  int mj2;
  double wjw[2];

  mj2 = 2 * mj;
  lj  = n / mj2;
  //cout << "n = " << n << endl;

# pragma omp parallel \
    shared ( a, b, c, d, lj, mj, mj2, sgn, w ) \
    private ( ambr, ambu, j, ja, jb, jc, jd, jw, k, wjw )

# pragma omp for nowait

  for ( j = 0; j < lj; j++ )
  {
    jw = j * mj;
    ja  = jw;
    jb  = ja;
    jc  = j * mj2;
    jd  = jc;
    //cout<<"mj = "<<mj<<endl;
    //cout << "jw,ja,jb = " << jw << endl;
    //cout << "jc,jd = " << jc << endl;
    wjw[0] = w[jw*2+0]; 
    wjw[1] = w[jw*2+1];
    //cout << "wjw[0] = " << wjw[0] << endl;
    //cout << "wjw[1] = " << wjw[1] << endl;

    if ( sgn < 0.0 ) 
    {
      wjw[1] = - wjw[1];
    }

    for ( k = 0; k < mj; k++ )
    {
      //cout << "in for k = " << k << endl;
      c[(jc+k)*2+0] = a[(ja+k)*2+0] + b[(jb+k)*2+0];
      c[(jc+k)*2+1] = a[(ja+k)*2+1] + b[(jb+k)*2+1];
      //cout << "c[" << (jc+k)*2+0 << "] = " << "a[" << (ja+k)*2+0 << "] + " << "b[" << (jb+k)*2+0 << "]" << endl;
     // cout << "c[" << (jc+k)*2+1 << "] = " << "a[" << (ja+k)*2+1 << "] + " << "b[" << (jb+k)*2+1 << "]" << endl;
      ambr = a[(ja+k)*2+0] - b[(jb+k)*2+0];
      ambu = a[(ja+k)*2+1] - b[(jb+k)*2+1];

     // cout << "ambr = " << "a[" << (ja+k)*2+0 << "] - " << "b[" << (jb+k)*2+0 << "]" << endl;
     // cout << "ambu = " << "a[" << (ja+k)*2+1 << "] - " << "b[" << (jb+k)*2+1 << "]" << endl;
      d[(jd+k)*2+0] = wjw[0] * ambr - wjw[1] * ambu;
      d[(jd+k)*2+1] = wjw[1] * ambr + wjw[0] * ambu;
      //cout << "d[" << (jd+k)*2+0 << "] = " << "wjw[0] * ambr - wjw[1] * ambu" << endl;
      //cout << "d[" << (jd+k)*2+1 << "] = " << "wjw[1] * ambr + wjw[0] * ambu" << endl;
    }
  }
  return;
}
void ccopy ( int n, double x[], double y[] )
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    y[i*2+0] = x[i*2+0];
    y[i*2+1] = x[i*2+1];
   }
  return;
}
void cfft2 ( int n, double x[], double y[], double w[], double sgn )
{
  int j;
  int m;
  int mj;
  int tgle;

   m = ( int ) ( log ( ( double ) n ) / log ( 1.99 ) );
   mj   = 1;
/*
  Toggling switch for work array.
*/
  tgle = 1;
  step ( n, mj, &x[0*2+0], &x[(n/2)*2+0], &y[0*2+0], &y[mj*2+0], w, sgn );

  if ( n == 2 )
  {
    return;
  }

  for ( j = 0; j < m - 2; j++ )
  {
    mj = mj * 2;
    if ( tgle )
    {
      step ( n, mj, &y[0*2+0], &y[(n/2)*2+0], &x[0*2+0], &x[mj*2+0], w, sgn );
      tgle = 0;
    }
    else
    {
      step ( n, mj, &x[0*2+0], &x[(n/2)*2+0], &y[0*2+0], &y[mj*2+0], w, sgn );
      tgle = 1;
    }
  }
/* 
  Last pass through data: move Y to X if needed.
*/
  if ( tgle ) 
  {
    ccopy ( n, y, x );
  }

  mj = n / 2;
  step ( n, mj, &x[0*2+0], &x[(n/2)*2+0], &y[0*2+0], &y[mj*2+0], w, sgn );

  return;
}
int main(int argc,char *argv[]) {
  //  freopen("test0.in","r",stdin);
  // freopen("test0b.out","w",stdout);
  float ti = 0;
  for(int c = 0; c < 10; c++)
  {
    wtime = omp_get_wtime();
    
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
    l = trans(n + m - 1);  // n次*m次不超过n+m-1次

    int *ans = new int[l];

    if(n == 1 && m == 1){
        cout << (s[0] - '0') * (t[0]-'0') << endl;
    }
    else{
    w = ( double * ) malloc (     l * sizeof ( double ) );
    x = ( double * ) malloc ( 2 * l * sizeof ( double ) );
    x1 = ( double * ) malloc ( 2 * l * sizeof ( double ) );
    y = ( double * ) malloc ( 2 * l * sizeof ( double ) );
    y2 = ( double * ) malloc ( 2 * l * sizeof ( double ) );
  long long int i;
  
# pragma omp parallel \
    shared ( l, x ) \
    private ( i)

# pragma omp for nowait
     for (i = 0; i < 2 * l; i = i + 2)
        {
          if(i/2 < n){
            x[i] = s[n-1-i/2] - '0';
          }
          else{
            x[i] = 0;
          }

          x[i+1] = 0;
        }

# pragma omp parallel \
    shared ( l, x1 ) \
    private ( i)

# pragma omp for nowait
        for (i = 0; i < 2 * l; i = i + 2)
        {
          if(i/2 < m){
            x1[i] = t[m-1-i/2] - '0';
          }
          else{
            x1[i] = 0;
          }

          x1[i+1] = 0;
        }
         /* for (int i = 0; i < 2 * l; i = i + 2 )
          {
          
            cout<<i<<" "<<x1[i]<<" "<<x1[i+1]<<endl;
            }cout<<endl;*/
        cffti ( l, w );
         sgn = + 1.0;
        cfft2 ( l, x, y, w, sgn );

       /* for (int i = 0; i < 2 * l; i = i + 2 )
          {
          
            cout<<i<<" "<<y[i]<<" "<<y[i+1]<<endl;
            }cout<<endl;*/
              cffti ( l, w );

         sgn = + 1.0;
        cfft2 ( l, x1, y2, w, sgn );

       /*  for (int i = 0; i < 2 * l; i = i + 2 )
          {
          
            cout<<i<<" "<<y2[i]<<" "<<y2[i+1]<<endl;
            }cout<<endl;*/
        #pragma omp parallel for
        for (int i = 0; i < 2 * l; i = i + 2 )
          {
            double temp = y[i];
            y[i] = y[i] * y2[i] - y[i+1]*y2[i+1];
            y[i+1] = temp*y2[i+1] + y[i+1]*y2[i];
        }

       /* for (int i = 0; i < 2 * l; i = i + 2 )
          {
          
            cout<<i<<" "<<y[i]<<" "<<y[i+1]<<endl;
            }cout<<endl;*/
        sgn = - 1.0;
        cfft2 ( l, y, x, w, sgn );

              for (int i = 0; i < 2 * l; i = i + 2 )
          {
              x[i] = x[i]/l;
              x[i+1] = x[i+1]/l;
          }

        /*
        for (int i = 0; i < 2 * l; i = i + 2 )
          {
            
            cout<<i<<" "<<x[i]<<" "<<x[i+1]<<endl;
            }cout<<endl;
        */
  
    #pragma omp parallel for
    for (int i = 0; i < l; i++) ans[i] = (int)(x[2*i] + 0.5);

    ans[l] = 0;  // error-prone :'l' -> '1'
    for (int i = 0; i < l; ++i) {
      ans[i + 1] += ans[i] / 10;
      ans[i] %= 10;
    }
    
    int p = l;
    
    for (; p && !ans[p]; --p)
      ;
    for (; ~p; putchar(ans[p--] + '0'))
      ;
        puts("");
    
    }

    wtime = omp_get_wtime() - wtime;
    cout << "time = " << wtime << endl;
    ti += (float)wtime;
     delete[] ans;
    delete[] s;
    delete[] t;
     free(w);
    free(x);
    free(x1);
    free(y);
    free(y2);
    }
    ti = ti/10.0;
    ofstream out;
    out.open("time.txt",ios::app);
    out<<ti<<endl;
    out.close();
   

    return 0;
}
