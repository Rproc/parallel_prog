#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <string.h>

#define _USE_MATH_DEFINES
#define SIZE 100000000

using namespace std;

int numDim = 2, n = 256;


struct data
{
	int tamanho;
	long double* vetDim;
}data;

typedef struct data Data;

struct dimensions
{
	Data* observations;
}dimensions;

typedef struct dimensions Observation;


int ReadFile1D (double*, int);
int ReadInputFile(double*, int);
double* KDE_Sequencial(double*, int, double);
void KDE_Mult(Observation*, double, long double*);
void createObservation(Observation*);



void createObservation(Observation* o){

	// int dimXY = 256;

	// X, Y
	for (int i = 0; i < numDim; ++i)
	{
		o->observations[i].vetDim = (long double*)malloc(sizeof(long double)*n);
		o->observations[i].tamanho = n;


	}

		for (int j = 0; j < n; ++j)
		{
			o->observations[0].vetDim[j] = j;
			o->observations[1].vetDim[j] = rand()%255;

		}

}

int main(int argc, char const *argv[])
{
 
	double bandwith = 0.01;

    srand( (unsigned)time(NULL) );

	Observation* o = (Observation*)malloc(sizeof(Observation));

	o->observations = (Data*)malloc(sizeof(Data)*numDim);

	createObservation(o);

	 // cout << "Tem que dar 0 " << o->observations[0].vetDim[0] << endl << "Tem que dar 1 "<< o->observations[1].vetDim[1] << endl;

	// cout << "Tem que dar alguma coisa " << o->observations[2].vetDim[0] << endl;


	long double* result = (long double*)malloc(n*sizeof(long double));

	for (int i = 0; i < n; ++i)
	{
		result[i] = 0.0;
	}

	
	KDE_Mult(o, bandwith, result);

	for (int i = 0; i < n; ++i)
	{
		cout << result[i] << " " ;

	}

	cout << endl;

	return 0;
}

// int ReadFile1D (double* input, int size){

//     std::fstream myfile ("/home/procopio/Documents/IC/Pro/file.txt", std::ios_base::in);
//     std::vector<double> numbers;

//     double x;
//     // char c;
    
//     while (myfile >> x){

//         numbers.push_back(x);
//         // numbers.push_back(y);
//     }

//     for(int i = 0; i < size && i < numbers.size(); i++){

//         input[i] = numbers[i];
//     }

//     // std::cout << "size: " << numbers.size() << endl;

//     return numbers.size();
// }


// int ReadInputFile(Observation* input, int size){

// 	std::fstream myfile("", std::ios_base::in);
// 	std::vector<double> numbers;

//     double x,y;
//     char c;
	
//     while (myfile >> x >> c >> y)
//     {
// 		numbers.push_back(x);
// 		numbers.push_back(y);
//     }

// 	for(int i = 0; i < size && i < numbers.size(); i++)
// 	{
// 		input[i] = numbers[i];
// 	}

// 	return numbers.size()/2;
// }


long double GaussianKernel(long double t)
{
    long double gaussian;

    gaussian = (1/sqrt(2 * M_PI)) * exp(-pow(t,2)/2);

    cout << "GAUSS  " <<gaussian << endl;

    return gaussian;
}

double* KDE_Sequencial(double* x, int size, double h){

	int i,j,k;
    double sum, prodKernel;
	double *pdf = new double[size];

	clock_t init, fim;
    double time_spent;
	init = clock();

    for(i = 0; i < size; i++)
    {
        sum = 0.0;
        for(j = 0; j < size; j++)
        {
            sum = sum + GaussianKernel( (x[i] - x[j] )/h) / h;
        
        }
        pdf[i] = sum/size;
		//printf("%f\n", pdf[i]);
    }

	fim = clock();
    time_spent = (double)(fim - init) / CLOCKS_PER_SEC;
	printf("KDE sequencial executado em %f segundos\n", time_spent);
	// printf("%f %f %f\n", pdf[0], pdf[size/2], pdf[size - 1]);


	return pdf;
}

void KDE_Mult(Observation* input, double h, long double* pdf){


	cout << "aaaa " << input->observations[0].tamanho << endl;

	for (int i = 0; i < n; ++i)
	{

		// cout << "[" << i << "]   Tam = " << input->observations[i].tamanho << endl;
		long double sum = 0.0;

		for (int j = 0; j < n; ++j)
		{

			long double prodKernel = 1.0;

			for (int k = 0; k < numDim; ++k)
			{
				cout << "Minus = " << (input->observations[k].vetDim[i] - input->observations[k].vetDim[j]) << endl;
				prodKernel = prodKernel * (GaussianKernel( (input->observations[k].vetDim[i] - input->observations[k].vetDim[j]) /h ) /h );
				// cout << "Prod [" <<i << ", " << j << "] = " << prodKernel << endl;
				// cout << "vetI = " << input->observations[k].vetDim[i] << " /// VetJ = " << input->observations[k].vetDim[j] << endl;
			}

			sum += prodKernel;
			 // cout << "Sum after Prod = " <<  sum << endl;
		}

		pdf[i] = sum/(double)n;
		// cout << pdf[i] << endl;
	}

}