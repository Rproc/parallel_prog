#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <string.h>

#define _USE_MATH_DEFINES
#define SIZE 100000000
#define DIMENSIONS 256

using namespace std;

int numDim = 3, t = 0, dimensionSize;


int ReadFile1D (double*, int);
int ReadInputFile(double*, int);
double* KDE_Sequencial(double*, int, double);
void KDE_Mult(double*, double, double*, int, int);
double* createInput(double*, int);


int main(int argc, char const *argv[]){

	// cout << "Num dim " << tam << endl;

	dimensionSize = atoi(argv[1]);

	t = atoi(argv[2]);

	// cout << "num threads = " << t << endl;

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


void KDE_Mult(double* input, double h, double* pdf, int inputSize, int size){


	FILE *f;
    f = fopen("result.txt", "a");


	clock_t ini, fim;
    double time_spent, time_spent_c;
	ini = clock();

	double init = omp_get_wtime();

	omp_set_num_threads(t);

	int i, j, k, id;
 	// double sum = 0.0, prodKernel = 1.0;
	#pragma omp parallel for private (i, j, k, id)
	for (i = 0; i < size; i++)
	{
		double sum = 0.0;

		for (j = 0; j < size; j++)
		{
			double prodKernel = 1.0;

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

	fprintf(f, "Tam = %d, num threads = %d\n", dimensionSize, t);
	fprintf(f, "KDE Open-MP executado em %f segundos\n", time_spent);

	fprintf(f, "KDE Open-MP executado em %f segundos (time C)\n", time_spent_c);

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
