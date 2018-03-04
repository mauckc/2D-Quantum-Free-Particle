//
//  freeparticle.c
//  This is a program intended to numerically integrate the schrodinger equation in 1-D
//  i*(dpsi/dx) = - (1/2)*( dpsi/dx)^2 + U(x)* psi(x)
//  where psi(x)= psireal(x) + i*psiimag(x)
//
//  This version implements a "split-step" Crank-Nicolson method.
//  We evolve our wave function in the position basis.
//  Then we fourier transform the wavefunction to evolve it in the momentum basis
//  Created by C. Ross Mauck on 1/23/15.
//
//  README
//  You must link this code with the fftw3 library. On Unix systems, link with -lfftw3 -lm.

#include "freeparticle.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <fftw3.h>
#include <complex>
#include <functional>
#include <ctime>

using namespace std;

#define N 256 // recommend 256 or above for stability
#define dt 0.01 // recommend 0.01 or below for stability
#define L 20.

#define RE 0
#define IM 1

#define t0 0.
#define tf 20.

#define div 20
#define desample 4 //desampling factor for output after calculation

#define AMPLITUDE 1.

double x = 0.0;
double y = 0.0;
double p = 0.0;
double px = 0.0;
double py = 0.0;
double kx = 0.0;
double ky = 0.0;

int num = 0; //for outfield iterator

double noutputs = (tf - t0)/(dt*div) ;

double t = 0.0; //this variable stores the time

double psi[2][2][N][N]; //this stores the wavefunction
double chi[2][2][N][N]; //this stores the wavefunction in fourier space

double dx = L/N;

double sigma = 0.73;
int slicenum = 0;
double A = 1.;

double realsum = 0.0;
double complexsum = 0.0;
double probabilitysum = 0.0;

void outputfield(int first)//outputs the field values
{
    static FILE *slicefield;
    static char name[500];

    sprintf(name,"./slices/slices_fields_%d.dat", first);
    slicefield=fopen(name,"w");

    double psiprob[N][N];

    for (int j = 0; j < N; j++)
    {
        for ( int i = 0 ; i < N; i++)
        {

        psiprob[j][i] = psi[0][RE][j][i] * psi[0][RE][j][i] + psi[0][IM][j][i] * psi[0][IM][j][i];

        }
    }

    for (int j = 0; j < N; j++)
    {
      for (int i = 0; i < N; i++)
      {
          if (i%desample==0){
          fprintf(slicefield,"%d  %d  %lf  %lf  %lf", i , j, psi[0][RE][j][i], psi[0][IM][j][i], psiprob[j][i]);
          fprintf(slicefield,"\n");
          }
        }
     }

    fclose(slicefield);
}

//void viewfield(int first)//shows fields
//{
  // Use openCV to draw a 2D image based on the relative propability of each pixel
//}

double potential(double x){
    double U = 50.0;
    return U * ( 1.0 - pow(cos(6*3.141592*(0.65*x-L/2)/L)/2,2.0));
}

void outputnorm(int first){

    double norm = 0.0;
    //norm +=

    return;
}

/*
void outputenergy(int first) //outputs energy value each slice.
{
    static FILE *sliceenergy;
    static char name[500];

    sprintf(name,"./slices/slices_energy_%d.dat", first);
    sliceenergy=fopen(name,"w");

    double energy = 0.0;

    for (int index = 0; index < N; index++)
    {
        if (index == 0)
        {
            energy +=(((-(psi[0][RE][index] - psi[0][IM][index])/2.) * (((psi[0][RE][N-1] + psi[0][IM][N-1] ) - 2 * (psi[0][RE][index] + psi[0][IM][index]) + (psi[0][RE][index+1] + psi[0][IM][index+1]))/(dx*dx)) + (psi[0][RE][index] - psi[0][IM][index]) * (psi[0][RE][index] + psi[0][IM][index]) * potential((index*dx) - (L/2)))/((psi[0][RE][index] - psi[0][IM][index]) * (psi[0][RE][index] + psi[0][IM][index])))*dx;
            //To account for normalization
            // (kinetic part) second derivative of wavefunction (I believe I will need to develop boundary conditions to deal with this
            //potential part
        }
        else if (index == N-1)
        {
            //congjugate wavefunction with -1/2 factor
            // (kinetic part) second derivative of wavefunction (I believe I will need to develop boundary conditions to deal with this
            //potential part
            energy += (((-(psi[0][RE][index] - psi[0][IM][index])/2.) * (((psi[0][RE][index-1] + psi[0][IM][index-1]) - 2 * (psi[0][RE][index] + psi[0][IM][index]) + (psi[0][RE][0] + psi[0][IM][0]))/(dx*dx)) + (psi[0][RE][index] - psi[0][IM][index]) * (psi[0][RE][index] + psi[0][IM][index]) * potential((index*dx) - (L/2)))/((psi[0][RE][index] - psi[0][IM][index]) * (psi[0][RE][index] + psi[0][IM][index])))*dx;

            //Last division to account for normalization
        }
        else
        {
            energy += (((-(psi[0][RE][index] - psi[0][IM][index])/2.) * (((psi[0][RE][index-1] + psi[0][IM][index-1]) - 2 * (psi[0][RE][index] + psi[0][IM][index]) + (psi[0][RE][index+1] + psi[0][IM][index+1]))/(dx*dx)) + (psi[0][RE][index] - psi[0][IM][index]) * (psi[0][RE][index] + psi[0][IM][index]) * potential((index*dx) - (L/2)))/((psi[0][RE][index] - psi[0][IM][index]) * (psi[0][RE][index] + psi[0][IM][index])))*dx;//To account for normalization
        }
    }


    fprintf(sliceenergy,"%lf", energy);
    fprintf(sliceenergy,"\n");

    fclose(sliceenergy);
}
*/

int main ()
{
    t=0.0;

    std::time_t result = std::time(nullptr);
    std::cout << std::asctime(std::localtime(&result))
              << result << " seconds since the Epoch\n";

    printf("number of output files: %lf", noutputs);

    fftw_complex *in, *out, *in2,*out2;
    fftw_plan plan, plan2;//plan will be foward and plan2 will be backward

    printf("\n***********************************\n\nSetting gaussian wavefunction... \n");

    for (int j = 0; j < N; j++){
        y = (j*dx) - (L/3.8);

        for (int i = 0; i < N; i++){

                x = (i*dx) - (L/3.8);
                kx = 20.0*3.141592/L;
                ky = 4.0*3.141592/L;
                A  = AMPLITUDE;


                psi[0][RE][j][i] = A*cos(kx*x+ky*y) * exp(-((x*x)/(4*sigma*sigma) + (y*y)/(4*sigma*sigma)));//Sets gaussian real part of psi
                psi[0][IM][j][i] = A*sin(kx*x+ky*y) * exp(-((x*x)/(4*sigma*sigma) + (y*y)/(4*sigma*sigma)));//Sets gaussian for the imaginary parts of psi
            }
    }

    printf("\nInitial conditions set \n");


        /*
    double realsum = 0.0;
    double complexsum = 0.0;
    double probabilitysum = 0.0;

    for ( int j = 0; j < N; j++ ){


    for (int index = 0; index < N; index++)
    {
        realsum += psi[0][RE][index];
        complexsum += psi[0][IM][index];
        probabilitysum += psi[0][RE][index] * psi[0][RE][index] + psi[0][IM][index] * psi[0][IM][index];
    }}

    printf("\n\n\nRealsum = %lf\nComplexsum = %lf\n\nProbabilitysum = %lf\n\n",realsum,complexsum,probabilitysum);*/
    outputfield(0);//init setting saved here

    //create both fftw plans to be used to time evolve our wavefunction
    //Forward DFT
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    plan = fftw_plan_dft_2d(N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    //Backward DFT
    in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    plan2 = fftw_plan_dft_2d(N, N, in2, out2, FFTW_BACKWARD, FFTW_ESTIMATE);

    printf("FFTw plans set!\n\n\n Ready to begin\n");

    printf("Beginning time-evolution\n");

    num=0;

    //time evolution begins
    while (t<tf)
    {

        for (int j = 0; j < N; j++){
            y = (j*dx - (L/2));
            for (int i = 0; i < N; i++){//update phase of position space wavefunction
            x = (i*dx - (L/2));

            psi[1][RE][j][i] = psi[0][RE][j][i] * cos(potential(x) * dt) + psi[0][IM][j][i]*sin(potential(x) * dt);
            psi[1][IM][j][i] = psi[0][IM][j][i] * cos(potential(x) * dt) - psi[0][RE][j][i]*sin(potential(x) * dt);

            }
        }


        for (int j = 0; j < N-1; j++){
            for (int i = 0; i < N; i++){  //load our FFTw array

            in[i+j*N][0] = psi[1][RE][j][i];
            in[i+j*N][1] = psi[1][IM][j][i];

            }
        }

        fftw_execute(plan);//transform now stored in out in FFTw format
        //this loop puts the transformed array in DFT output


            for (int j = 0; j < N-1; j++){
                for (int i = 0; i < N; i++){

                    chi[0][RE][j][i] = out[i+j*N][0];
                    chi[0][IM][j][i] = out[i+j*N][1];
                }
            }


        int index = 0;
        for (int j = 0; j < N-1; j++)
        {

            py = ((2*3.145926535)/L) * (( (j + (N/2)) % N) - N/2);

            for (int i = 0; i < N; i++)//here we update the phases in momentum space
            {
                px = ((2*3.145926535)/L) * (( (i + (N/2)) % N) - N/2);

                chi[1][RE][j][i] = chi[0][IM][j][i]*sin((dt*(px*px+py*py))/2) + chi[0][RE][j][i]*cos((dt*(px*px+py*py))/2);
                chi[1][IM][j][i] = chi[0][IM][j][i]*cos((dt*(px*px+py*py))/2) - chi[0][RE][j][i]*sin((dt*(px*px+py*py))/2);
            }

        }

        for (int j = 0; j < N-1; j++){
            for (int i = 0; i < N; i++){
                //load our FFTw array
            in2[i+j*N][0] = chi[1][RE][j][i];
            in2[i+j*N][1] = chi[1][IM][j][i];
            }
        }
        fftw_execute(plan2);

            for (int j = 0; j < N-1; j++){
                for (int i = 0; i < N; i++){
                        //this loop puts the transformed array in DFT output

            psi[0][RE][j][i] = out2[i+j*N][0];
            psi[0][IM][j][i] = out2[i+j*N][1];
                }
            }
        //this loop accounts for unnormalized DFT after fwd and bkwd trnsfms
            for (int j = 0; j < N; j++){
                for (int i = 0; i < N; i++){  //load our FFTw array

            psi[0][RE][j][i] = psi[0][RE][j][i]/(N*N);
            psi[0][IM][j][i] = psi[0][IM][j][i]/(N*N);
                }
            }


        t += dt;
        num++;
        printf("*** program time: %lf \n",t);
        if (num % div == 0)
        {
            outputfield(slicenum);
            //outputenergy(slicenum);
            slicenum++;

        }


         }

    printf("*** time-evolution has been terminated ***\n");
    fftw_destroy_plan(plan);
    fftw_destroy_plan(plan2);

    fftw_free(in); fftw_free(in2); fftw_free(out); fftw_free(out2);

    std:time_t resultold = result;
    result = std::time(nullptr);
    std::cout << std::asctime(std::localtime(&result))
              << result - resultold << " seconds since the Epoch\n";

    printf("\n____ number of output files: %lf ____", noutputs);


    return 0;



}
