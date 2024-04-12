#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <string.h>

#define _USE_MATH_DEFINES
#define SIZE 100000000
#define DIMENSIONS 64

using namespace std;

int numDim = 3;


int ReadFile1D (double*, int);
int ReadInputFile(double*, int);
double* KDE_Sequencial(double*, int, double);
void KDE_Mult(double*, double, double*, int, int);
double* createInput(double*, int);


int main(int argc, char const *argv[]){

	// cout << "Num dim " << tam << endl;

	int dimensionSize = atoi(argv[1]);


	double bandwith = 0.5;
	srand( (unsigned)time(NULL) );

	int size = dimensionSize*dimensionSize;
	int inputSize = size * numDim;

	// cout << "size " << size << "\t\t" << "inputSize = " << inputSize << endl;

	double* input = (double*)malloc(inputSize * sizeof(double));
	double* pdf = (double*)malloc(size * sizeof(double));

	memset(input, 0, inputSize * sizeof(double));
	memset(pdf, 0, size * sizeof(double));

	input = createInput(input, inputSize);

	// cout << "Input pos[0] = " << input[3] << endl;

	KDE_Mult(input, bandwith, pdf, inputSize, size);

	// for (int i = 0; i < size; ++i)
	// {
	// 	cout << pdf[i] << " " ;
	// }

	cout << endl;

	return 0;
}

double* createInput(double* input, int inputSize){

	for (int i = 0; i < inputSize; i++) {

		input[i] = rand()%256;
	}

	return input;
}

long double GaussianKernel(long double t)
{
	long double gaussian;

	gaussian = (1/sqrt(2 * M_PI)) * exp(-pow(t,2)/2);

	// cout << "GAUSS  " <<gaussian << endl;

	return gaussian;
}

// double* KDE_Sequencial(double* x, int size, double h){
//
// 	int i,j,k;
// 	double sum, prodKernel;
// 	double *pdf = new double[size];
//
// 	clock_t init, fim;
// 	double time_spent;
// 	init = clock();
//
// 	for(i = 0; i < size; i++)
// 	{
// 		sum = 0.0;
// 		for(j = 0; j < size; j++)
// 		{
// 			sum = sum + GaussianKernel( (x[i] - x[j] )/h) / h;
//
// 		}
// 		pdf[i] = sum/size;
// 		//printf("%f\n", pdf[i]);
// 	}
//
// 	fim = clock();
// 	time_spent = (double)(fim - init) / CLOCKS_PER_SEC;
// 	printf("KDE sequencial executado em %f segundos\n", time_spent);
// 	// printf("%f %f %f\n", pdf[0], pdf[size/2], pdf[size - 1]);
//
//
// 	return pdf;
// }

void KDE_Mult(double* input, double h, double* pdf, int inputSize, int size){

	clock_t ini, fim;
    double time_spent, time_spent_c;
	ini = clock();

	double init = omp_get_wtime();

	int i, j, k;
 	double sum = 0.0, prodKernel = 1.0;
	#pragma omp parallel for private (i, j, k, sum, prodKernel)
	for (i = 0; i < size; i++)
	{
		 sum = 0.0;

		// cout << "HEY" << endl;
		for (j = 0; j < size; j++)
		{
			 prodKernel = 1.0;
			// cout << "HEY" << endl;

			for (k = 0; k < numDim; k++)
			{
				#pragma omp atomic
				prodKernel *= ( GaussianKernel( (input[k * size + i] - input[k * size + j]) /h) /h );
			}
			#pragma omp atomic
			sum += prodKernel;
			// cout << "Sum after Prod = " <<  sum << endl;
		}

		pdf[i] = sum/(double)size;
	}
	double end = omp_get_wtime();

	fim = clock();
    time_spent = (double)(end - init);
	time_spent_c = (double)(fim - ini) / CLOCKS_PER_SEC;

	printf("KDE Open-MP executado em %f segundos\n", time_spent);

	printf("KDE Open-MP executado em %f segundos (time C)\n", time_spent_c);

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
//
// 	std::fstream myfile("", std::ios_base::in);
// 	std::vector<double> numbers;
//
//     double x,y;
//     char c;
//
//     while (myfile >> x >> c >> y)
//     {
// 		numbers.push_back(x);
// 		numbers.push_back(y);
//     }
//
// 	for(int i = 0; i < size && i < numbers.size(); i++)
// 	{
// 		input[i] = numbers[i];
// 	}
//
// 	return numbers.size()/2;
// }
//
