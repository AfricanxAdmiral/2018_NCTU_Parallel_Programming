/**********************************************************************
 * DESCRIPTION:
 *   Serial Concurrent Wave Equation - C Version
 *   This program implements the concurrent wave equation
 *********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAXPOINTS 1000000
#define MAXSTEPS 1000000
#define MINPOINTS 20
#define PI 3.14159265
#define THREAD 216

void check_param(void);
__global__ void update(int tpoints, int nsteps, float *values);
void printfinal(void);

int nsteps,                     /* number of time steps */
    tpoints,                    /* total points along string */
    rcode;                      /* generic return code */
float  values[MAXPOINTS+2];     /* values at time t */
       //oldval[MAXPOINTS+2],   /* values at time (t-dt) */
       //newval[MAXPOINTS+2];   /* values at time (t+dt) */

/**********************************************************************
 *      Checks input values from parameters
 *********************************************************************/
void check_param(void)
{
   char tchar[20];

   /* check number of points, number of iterations */
   while ((tpoints < MINPOINTS) || (tpoints > MAXPOINTS)) {
      printf("Enter number of points along vibrating string [%d-%d]: "
           ,MINPOINTS, MAXPOINTS);
      scanf("%s", tchar);
      tpoints = atoi(tchar);
      if ((tpoints < MINPOINTS) || (tpoints > MAXPOINTS))
         printf("Invalid. Please enter value between %d and %d\n",
                 MINPOINTS, MAXPOINTS);
   }
   while ((nsteps < 1) || (nsteps > MAXSTEPS)) {
      printf("Enter number of time steps [1-%d]: ", MAXSTEPS);
      scanf("%s", tchar);
      nsteps = atoi(tchar);
      if ((nsteps < 1) || (nsteps > MAXSTEPS))
         printf("Invalid. Please enter value between 1 and %d\n", MAXSTEPS);
   }

   printf("Using points = %d, steps = %d\n", tpoints, nsteps);

}
/**********************************************************************
 *     Update all values along line a specified number of times
 *********************************************************************/
__global__ void update(int tpoints, int nsteps, float *values)
{
        int i;
        int j = blockIdx.x * blockDim.x + threadIdx.x;
        float x, fac, tmp;
        __shared__ float value, oldval, newval;
        /* Calculate initial values based on sine curve */
        fac = 2.0 * PI;
        tmp = tpoints - 1;
        x = (float)(j-1)/tmp;
        value = sin (fac * x);
        oldval = value;
        //float dtime, c, dx, tau;
        float sqtau;
        //dtime = 0.3;
        //c = 1.0;
        //dx = 1.0;
        //tau = (c * dtime / dx);
        sqtau = 0.3 * 0.3;
        /* Update values for each time step */
        if ((j == 1) || (j == tpoints))
                value = 0.0;
        else {
                for (i = 1; i<= nsteps; i++) {
                        //newval = 1.82 * value - oldval;
                        newval = (2.0 * (1.0 - sqtau) * value) - oldval;
                        oldval = value;
                        value = newval;
                }
        }
        values[j] = value;
}
/**********************************************************************
 *     Print final results
 *********************************************************************/
void printfinal()
{
   int i;

   for (i = 1; i <= tpoints; i++) {
      printf("%6.4f ", values[i]);
      if (i%10 == 0)
         printf("\n");
   }
}
/**********************************************************************
 *      Main program
 *********************************************************************/
int main(int argc, char *argv[])
{
        sscanf(argv[1],"%d",&tpoints);
        sscanf(argv[2],"%d",&nsteps);
        check_param();

        float *cuda_values;
        int size = (tpoints+1) * sizeof(float);
        cudaMalloc(&cuda_values, size);

        printf("Initializing points on the line...\n");
        printf("Updating all points for all time steps...\n");
        update<<<(tpoints+THREAD)/THREAD, THREAD>>>(tpoints, nsteps, cuda_values);
        cudaMemcpy(values, cuda_values, size, cudaMemcpyDeviceToHost);
        printf("Printing final results...\n");
        printfinal();
        printf("\nDone.\n\n");

        cudaFree(cuda_values);

        return 0;
}
