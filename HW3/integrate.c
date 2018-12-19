#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define PI 3.1415926535

int main(int argc, char **argv) 
{
  int rank;
  int size;
  int source;
  int dest = 0;
  int tag = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  long long i, num_intervals;
  double rect_width, area, sum, my_sum, x_middle; 

  sscanf(argv[1],"%llu",&num_intervals);

  rect_width = PI / num_intervals;

  my_sum = 0;
  for(i = rank + 1; i < num_intervals + 1; i += size) {

    /* find the middle of the interval on the X-axis. */ 

    x_middle = (i - 0.5) * rect_width;
    area = sin(x_middle) * rect_width; 
    my_sum += area;
  } 

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Reduce(&my_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, dest, MPI_COMM_WORLD);

  if (!rank)
    printf("The total area is: %f\n", (float)sum);

  MPI_Finalize();

  return 0;
}   
