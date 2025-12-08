#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "common.h"
#include <inttypes.h>
#include <time.h>


#define PI 3.1415926535


void usage(int argc, char** argv);
double calcPi_Serial(int num_steps);
double calcPi_P1(int num_steps);
double calcPi_P2(int num_steps);

double half_circle(double x) 
{
    // Height of half unit circle at point x with origin at 0
    return sqrt(1 - pow(x, 2));
}

int main(int argc, char** argv)
{
    // get input values
    uint32_t num_steps = 100000;
    if(argc > 1) {
        num_steps = atoi(argv[1]);
    } else {
        usage(argc, argv);
        printf("using %"PRIu32"\n", num_steps);
    }
    fprintf(stdout, "The first 10 digits of Pi are %0.10f\n", PI);


    // set up timer
    uint64_t start_t;
    uint64_t end_t;
    InitTSC();


    // calculate in serial 
    start_t = ReadTSC();
    double Pi0 = calcPi_Serial(num_steps);
    end_t = ReadTSC();
    printf("Time to calculate Pi serially with %"PRIu32" steps is: %g\n",
           num_steps, ElapsedTime(end_t - start_t));
    printf("Pi is %0.10f\n", Pi0);
    
    // calculate in parallel with integration
    start_t = ReadTSC();
    double Pi1 = calcPi_P1(num_steps);
    end_t = ReadTSC();
    printf("Time to calculate Pi in // with %"PRIu32" steps is: %g\n",
           num_steps, ElapsedTime(end_t - start_t));
    printf("Pi is %0.10f\n", Pi1);


    // calculate in parallel with Monte Carlo
    start_t = ReadTSC();
    double Pi2 = calcPi_P2(num_steps);
    end_t = ReadTSC();
    printf("Time to calculate Pi in // with %"PRIu32" guesses is: %g\n",
           num_steps, ElapsedTime(end_t - start_t));
    printf("Pi is %0.10f\n", Pi2);

    
    return 0;
}


void usage(int argc, char** argv)
{
    fprintf(stdout, "usage: %s <# steps>\n", argv[0]);
}

double calcPi_Serial(int num_steps)
{
    double pi = 0.0;
    double dx = 2.0 / num_steps;

    for (int i = 0; i < num_steps; i++) {
        pi += half_circle(-1.0 + (dx * i)) * dx;
    }
    return 2 * pi;
}

double calcPi_P1(int num_steps)
{
    double pi = 0.0;
    double dx = 2.0 / num_steps;

    #pragma omp parallel for
    for (int i = 0; i < num_steps; i++) {
        #pragma omp automic
        pi += half_circle(-1 + (dx * i)) * dx;
    }
    return 2 * pi;
}

double calcPi_P2(int num_steps)
{
    double pi = 0.0;
    double count = 0.0;

    #pragma omp parallel for
    for (unsigned int i = 0; i < num_steps; i++) {
        double x = (2 * ((double)(rand_r(&i)) / (double)(RAND_MAX))) - 1.0;
        double y = (2 * ((double)(rand_r(&i)) / (double)(RAND_MAX))) - 1.0;
        
        if (pow(x, 2) + pow(y, 2) <= 1.0) {
            #pragma omp atomic
            count += 1;
        }
    }
    
    pi = ((double)(count * 4)) / ((double)num_steps);
    return pi;
}
