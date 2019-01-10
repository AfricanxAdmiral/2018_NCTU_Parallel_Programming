#include <cstdio>
#include <cmath>
#include <complex>
#include <cstring>
#include<iostream>
#include <pthread.h>
#include <omp.h>
#include<fstream>
#include<cstdlib>
using namespace std;
const double PI(acos(-1.0));
typedef complex<double> C;

const int N = (1 << 30);

int thread_count = 4; 
pthread_t *thread_handles;
char *s;
char *t;
  int n, m, l;
  int mj;
  double *a,*b,*c,*d;
  int *ans;
  void* step_thread(void* rank);
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
/*
# pragma omp parallel \
    shared ( aw, n, w ) \
    private ( arg, i )

# pragma omp for nowait
*/
  for ( i = 0; i < n2; i++ )
  {
    arg = aw * ( ( double ) i );
    w[i*2+0] = cos ( arg );
    w[i*2+1] = sin ( arg );
  }
  return;
}

void* step_thread(void* rank){
  double ambr;
  double ambu;
  long j;
  int ja;
  int jb;
  int jc;
  int jd;
  int jw;
  int k;
  int lj;
  int mj2;
  double wjw[2];
  long my_rank = (long )rank;

  mj2 = 2 * mj;
  lj  = l / mj2;

  long start,end;
  long workload = lj/thread_count;
  if(my_rank == 0){
    start = 0;
    end = start + workload;
  }
  else if(my_rank == thread_count - 1){
    start = my_rank * workload;
    end = lj;
  }
  else{
    start = my_rank * workload;
    end = start + workload;
  }

  for ( j = start; j < end; j++ )
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
  return NULL;
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
  int tgle;

   m = ( int ) ( log ( ( double ) n ) / log ( 1.99 ) );
   mj   = 1;
/*
  Toggling switch for work array.
*/
  tgle = 1;
  a = &x[0*2+0];
  b = &x[(n/2)*2+0];
  c = &y[0*2+0];
  d = &y[mj*2+0];

  for(long i = 0; i < thread_count; i++){
    pthread_create(&thread_handles[i],NULL,step_thread,(void*)i);
  }

  for(long i = 0; i < thread_count; i++){
    pthread_join(thread_handles[i],NULL);
  }

  if ( n == 2 )
  {
    return;
  }

  for ( j = 0; j < m - 2; j++ )
  {
    mj = mj * 2;
    if ( tgle )
    {
        a = &y[0*2+0];
        b = &y[(n/2)*2+0];
        c = &x[0*2+0];
        d = &x[mj*2+0];

      for(long i = 0; i < thread_count; i++){
        pthread_create(&thread_handles[i],NULL,step_thread,(void*)i);
      }

      for(long i = 0; i < thread_count; i++){
        pthread_join(thread_handles[i],NULL);
      }
      tgle = 0;
    }
    else
    {
        a = &x[0*2+0];
        b = &x[(n/2)*2+0];
        c = &y[0*2+0];
        d = &y[mj*2+0];
        for(long i = 0; i < thread_count; i++){
          pthread_create(&thread_handles[i],NULL,step_thread,(void*)i);
        }

        for(long i = 0; i < thread_count; i++){
          pthread_join(thread_handles[i],NULL);
        }
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
    a = &x[0*2+0];
        b = &x[(n/2)*2+0];
        c = &y[0*2+0];
        d = &y[mj*2+0];
        for(long i = 0; i < thread_count; i++){
          pthread_create(&thread_handles[i],NULL,step_thread,(void*)i);
        }

        for(long i = 0; i < thread_count; i++){
          pthread_join(thread_handles[i],NULL);
        }

  return;
}

void* initialize(void* rank){
	/*  
# pragma omp parallel \
    shared ( l, x ) \
    private ( i)

# pragma omp for nowait
*/
	
	long my_rank = (long)rank;
	//cout << "my_rank = " << my_rank << endl;
	//cout << "rank finish\n";
	long start,end;
  long workload = l/thread_count;
  //cout << "workload = " << workload << endl;
  if(my_rank == 0){
    start = 0;
    end = start + workload;
  }
  else if(my_rank == thread_count - 1){
    start = 2 * my_rank * workload;
    end = l;
  }
  else{
    start = my_rank * workload;
    end = start + workload;
    start = 2 * start;
  }
  //cout << "l = " << l << endl;
  /*
  if(my_rank == 3){
  	cout << "start = " << start << endl;
  	cout << "end = " << end << endl;
  }
  */
  //cout << "start execute\n";
    for (int i = start; i < 2 * end; i = i + 2)
    {
    	//cout << "in n i = " << i << endl;
        if(i/2 < n){
            x[i] = s[n-1-i/2] - '0';
        }
        else{
            x[i] = 0;
        }
        x[i+1] = 0;
    }
/*
# pragma omp parallel \
    shared ( l, x1 ) \
    private ( i)

# pragma omp for nowait
*/
    for (int i = start; i < 2 * end; i = i + 2)
    {
    	//cout << "in m i = " << i << endl;
        if(i/2 < m){
            x1[i] = t[m-1-i/2] - '0';
        }
        else{
            x1[i] = 0;
        }
        x1[i+1] = 0;
    }
}

void* multiply(void* rank){
	/*  
# pragma omp parallel \
    shared ( l, x ) \
    private ( i)

# pragma omp for nowait
*/
	
	long my_rank = (long)rank;
	//cout << "my_rank = " << my_rank << endl;
	//cout << "rank finish\n";
	long start,end;
  long workload = l/thread_count;
  //cout << "workload = " << workload << endl;
  if(my_rank == 0){
    start = 0;
    end = start + workload;
  }
  else if(my_rank == thread_count - 1){
    start = 2 * my_rank * workload;
    end = l;
  }
  else{
    start = my_rank * workload;
    end = start + workload;
    start = 2 * start;
  }
  
  for (int i = start; i < 2 * end; i = i + 2 )
    {
        double temp = y[i];
        y[i] = y[i] * y2[i] - y[i+1]*y2[i+1];
        y[i+1] = temp*y2[i+1] + y[i+1]*y2[i];
    }
}

void* final(void* rank){
	/*  
# pragma omp parallel \
    shared ( l, x ) \
    private ( i)

# pragma omp for nowait
*/
	
	long my_rank = (long)rank;
	//cout << "my_rank = " << my_rank << endl;
	//cout << "rank finish\n";
	long start,end,start1,end1;
  long workload = l/thread_count;
  //cout << "workload = " << workload << endl;
  if(my_rank == 0){
    start = 0;
    end = start + workload;
    start1 = 0;
    end1 = start1 + workload;
  }
  else if(my_rank == thread_count - 1){
    start = 2 * my_rank * workload;
    end = l;
    start1 = my_rank * workload;
    end1 = l;
  }
  else{
    start = my_rank * workload;
    end = start + workload;
    start = 2 * start;
    start1 = my_rank * workload;
    end1 = start1 + workload;
  }
  
  for (int i = start; i < 2 * end; i = i + 2 )
    {
        x[i] = x[i]/l;
        x[i+1] = x[i+1]/l;
    }
    
    //#pragma omp parallel for
    for (int i = start1; i < end1; i++) ans[i] = (int)(x[2*i] + 0.5);
}

int main(int argc,char *argv[]) {
  //  freopen("test0.in","r",stdin);
  // freopen("test0b.out","w",stdout);
    thread_count = atoi(argv[3]);
    wtime = omp_get_wtime ( );
    thread_handles = (pthread_t*)malloc(thread_count * sizeof(pthread_t));

    s = new char[N];
    t = new char[N];


  
    FILE* fin = fopen(argv[1],"r");
    fscanf(fin,"%s",s);
    fclose(fin);
	
    fin = fopen(argv[2],"r");
    fscanf(fin,"%s",t);
    fclose(fin);
	
    
    n = strlen(s);
    m = strlen(t);
    l = trans(n + m - 1);  // n次*m次不超过n+m-1次

    ans = new int[l];
	
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
	for(long i = 0; i < thread_count; i++){
		//cout << "start create thread " << i << endl;
   	 	pthread_create(&thread_handles[i],NULL,initialize,(void*)i);
   	 	//cout << "create thread " << i << " finish" << endl;
  	}

	for(long i = 0; i < thread_count; i++){
    pthread_join(thread_handles[i],NULL);
  	}
  	
  	
    cffti ( l, w );

    sgn = + 1.0;
    cfft2 ( l, x, y, w, sgn );

    sgn = + 1.0;
    cfft2 ( l, x1, y2, w, sgn );
    
    for(long i = 0; i < thread_count; i++){	
   	 	pthread_create(&thread_handles[i],NULL,multiply,(void*)i);

  	}

	for(long i = 0; i < thread_count; i++){
    pthread_join(thread_handles[i],NULL);
  	}

    sgn = - 1.0;
    cfft2 ( l, y, x, w, sgn );
	
	for(long i = 0; i < thread_count; i++){	
   	 	pthread_create(&thread_handles[i],NULL,final,(void*)i);

  	}

	for(long i = 0; i < thread_count; i++){
    pthread_join(thread_handles[i],NULL);
  	}

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


    wtime = omp_get_wtime ( ) - wtime;
    cout << "time = " << wtime << endl;

    delete[] ans;
    delete[] s;
    delete[] t;
     free(w);
    free(x);
    free(x1);
    free(y);
    free(y2);
  
    ofstream out;
    out.open("time_pthread.txt",ios::app);
    out<<wtime<<" ";
    out.close();

    return 0;
}
