#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int my_pc;
int my_done;
int size;
long long int limit;

void prime(int rank) {
  int n;
  my_pc = 0;
  for (n=11+2*rank; n<=limit; n=n+2*size) {
    int i,squareroot, flag = 1;
    if (n>10) {
      squareroot = (int) sqrt(n);
      for (i=3; i<=squareroot; i=i+2)
        if ((n%i)==0){
          flag = 0;
          break;
        }
      if(flag){
        my_pc++;
        my_done = n;
      }
    }
  }
}

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  int rank;
  int source;
  int dest = 0;
  int tag = 0;
  int messeage;
  MPI_Status status;
  MPI_Request request;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  int pc;
  int foundone; /* most recent prime found */
  long long int batch, remain;

  sscanf(argv[1],"%llu",&limit);	
  if (!rank)
    printf("Starting. Numbers to be scanned= %lld\n",limit);

  batch = limit / size;
  pc=4;     /* Assume (2,3,5,7) are counted here */

  prime(rank);
  int a[2] = {my_pc, my_done};

  /*
  for (n=11+2*rank; n<=limit; n=n+2*size) {
    if (isprime(n)) {
      my_pc++;
      my_done = n;
    }			
  }
  */

  /*
  if (!rank){
    foundone = my_done;
    pc += my_pc;
    for (source = 1; source < size; source++){
      MPI_Recv(&a, 2, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
      //MPI_Wait(&request, &status);
      foundone = foundone < a[1] ? a[1] : foundone;
      pc += a[0];
    }
  }
  else{
    MPI_Send(&a, 2, MPI_INT, dest, tag, MPI_COMM_WORLD);
  }
  */

  MPI_Barrier(MPI_COMM_WORLD);

  /* 
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  printf("Host : %s\n", processor_name);
  */

  MPI_Reduce(&my_done, &foundone, 1, MPI_INT, MPI_MAX, dest, MPI_COMM_WORLD);
  MPI_Reduce(&my_pc, &pc, 1, MPI_INT, MPI_SUM, dest, MPI_COMM_WORLD);

  if (!rank)
    printf("Done. Largest prime is %d Total primes %d\n",foundone,pc + 4);
  MPI_Finalize();
  return 0;
}
