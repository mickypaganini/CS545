// --------------------------------------------------------------------//
// Michela Paganini                      							   //
// CS 545 - Introduction to Data Mining 							   // 
// Yale University                 								       //
// October 6, 2014                  							       //
// Homework 2                            							   //
// "MINIMIZATION VIA STEEPEST DESCENT"                                 //
// GOAL: Solve Ax=y by minimizing f(x)=||Ax-y||^2 by steepest descent  //
// --------------------------------------------------------------------//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// ---------------------------------------//
//                Prototypes              //
// ---------------------------------------//
void dumb_solve(double *a, double *y, int n, double eps, int numit, double *x, int *niter, double *discreps);
void print_arr(FILE *f, double * arr, int n);
void print_mat(FILE *f, double * mat, int n_row, int n_col);
void matrix_multiply(double *a, double *b, double *ab, int rows_a, int cols_a, int rows_b, int cols_b);
double length_sq(double *arr, int n);
void grad(double *DelF, double *a, double *ax_y, int n);
void transpose(double *mat, double *trans, int n_row, int n_col);
double line_search(double *a, double *DelF, int n);
// ---------------------------------------//
//                 Functions              //
// ---------------------------------------//
void dumb_solve(double *a, double *y, int n, double eps, int numit, double *x, int *niter, double *discreps){
	// ------------------------------------------------------------------------------------//
	// DEFINITIONS:
	// Inputs:
	// a = nxn matrix of the system, given as real nxn array
	// y = real array of size n representing the transformed vector
	// n = dimensions, size of system
	// eps = relative accuracy to which system is to be solved
	// numit = max number of iterations to be performed
	// Outputs:
	// x = real array of size n representing solution to the system
	// niter = number of iterations actually performed
	// discreps = array of length niter representing discrepancies after each iteration
	// ------------------------------------------------------------------------------------//
	
	while (*niter < numit){ // break once niter = numit 

		// Find vector (Ax-y):
		double *ax = malloc(n *sizeof(double));
		double *ax_y = malloc(n *sizeof(double));
		matrix_multiply(a, x, ax, n, n, n, 1); // calculates a*x at initial point x, puts it in ax

		int i;
		for (i = 0; i < n; ++i){ 
			ax_y[i] = ax[i] - y[i]; // calculates ax-y at initial point x
		}
		
		// Find ∇F(x) = 2(A)^T * (Ax-y) at initial point x:
		double *DelF = calloc(n, sizeof(double));
		grad(DelF, a, ax_y, n); // calls function that calculates ∇F(x) = 2(A)^T * (Ax-y), stores it into DelF
		
		// Find step: 
		double step = -line_search(a, DelF, n); // using line search method

		// Find the next position x:
		int j;
		for (j = 0; j < n; ++j){
			x[j] += step * DelF[j]; // x_{n+1} = x_n + step * ∇F(x)
		}
		
		++(*niter); // done with the iteration, so increase the count

		// Find the discrepancies:
		matrix_multiply(a, x, ax, n, n, n, 1); // calculates ax at new point x
		int k;
		for (k = 0; k < n; ++k){ 
			ax_y[k] = ax[k] - y[k]; // calculates ax-y at new point x
		}
		double F = length_sq(ax_y, n); // calculates F(x) = ||Ax-y||^2 at the new position
		discreps[*niter-1] = F; // stores ||Ax-y||^2 as an element of the array discreps (-1 because array index starts at 0)

		// If relative accuracy is already reached, stop the iteration:
		if (F < eps){
			break;  
		}
		
		// Free allocated memory:		
		free(ax); free(ax_y); free(DelF);
		
	}	
	
}
// ---------------------------------------//
void print_arr(FILE *f, double * arr, int n){ 
// function that prints an array of size n 
	int i;
	printf("[ ");
	fprintf(f, "[ ");

	for (i = 0; i < (n - 1); ++i){
		printf("%g , ", arr[i]);
		fprintf(f, "%g , ", arr[i]);
	}
	printf("%g ]\n", arr[i]);
	fprintf(f, "%g ]\n", arr[i]);
}
// ---------------------------------------//
void print_mat(FILE *f, double * mat, int n_row, int n_col){
// function that prints a matrix of size n_row x n_col
	int i, j;
	for (i = 0; i < n_row; ++i){
		for (j = 0; j < n_col; ++j){
			printf("%g  ", mat[i * n_col + j]);
			fprintf(f, "%g  ", mat[i * n_col + j]);
		}
		printf("\n");
		fprintf(f, "\n");
	}
	printf("\n");
	fprintf(f, "\n");
}
// ---------------------------------------//
void matrix_multiply(double *a, double *b, double *ab, int rows_a, int cols_a, int rows_b, int cols_b){
// function that multiplies matrices a and b, storing product in matrix ab
	if (cols_a != rows_b){ 
		printf("ERROR! The number of columns of A must equal the number of rows of B.");
	}
	else{
		int row, col, index, i;
		// Initialize matrix ab to have all entries = 0, so that values can be incremented by loop:
		for (i = 0; i < (rows_a * cols_b); ++i){ 
			ab[i] = 0.0;
		}
		// Find entries of ab:
		for (row = 0; row < rows_a; ++row){
			for (col = 0; col < cols_b; ++col){
				for (index = 0; index < cols_a; ++index){
					ab[row * cols_b + col] += a[row * cols_a + index] * b[index * cols_b + col];
				}
			}
		}
	}
}
// ---------------------------------------//
double length_sq(double *arr, int n){
// function that calculates the norm squared of an array
	int i;
	double sum_sq = 0;
	for (i = 0; i < n; ++i){
		sum_sq += pow(arr[i], 2.0);
	}
	return sum_sq;
}
// ---------------------------------------//
void grad(double *DelF, double *a, double *ax_y, int n){
// function that calculates ∇F(x) = 2(A)^T * (Ax-y)
	int i, j;
	double *AT = calloc(n * n, sizeof(double));
	double *twoAT = calloc(n * n, sizeof(double));
	transpose(a, AT, n, n); // calculates (A)^T to update AT

	// Find 2*(A)^T:
	for (i = 0; i < n; ++i){
		for (j = 0; j < n; ++j){
			twoAT[i * n + j] = 2 * AT[i * n + j];
		}
	}

	// Find 2*(A)^T * (Ax-y):
	matrix_multiply(twoAT, ax_y, DelF, n, n, n, 1); // stores it into DelF
	// Free allocated memory:
	free(AT); free(twoAT);
}	
// ---------------------------------------//
void transpose(double *mat, double *trans, int n_row, int n_col){
// function that calculates the transpose of a matrix
	int i,j;
	for (i = 0; i < n_col; ++i){
		for (j = 0; j < n_row; ++j){
			trans[i * n_row + j] = mat[j * n_col + i];
		}
	}
}
// ---------------------------------------//
double line_search(double *a, double *DelF, int n){
// function that calculates step via formula: alpha = ||∇F(x)||^2 / [2 * ||A*∇F(x)||^2]
	double *aDelF = calloc(n, sizeof(double));
	matrix_multiply(a, DelF, aDelF, n, n, n, 1); // computes A*∇F(x) and stores it in aDelF
	double alpha = length_sq(DelF, n) / (2 * length_sq(aDelF, n));
	free(aDelF);
	return alpha;
}
// --------------------------------------------//
//                  Main                       //
// 				  (Check)					   //
// --------------------------------------------// 
int main(){
	FILE *f;
	f = fopen("hw2.txt", "w"); // output file

	// Initialize variables:
	int N[4] = {6, 3, 4, 10}; // dimensions of diagonal matrix A (nxn); four different experiments
	int numit = 1000; // max number of iterations = 1000
	double eps = pow(10.0, -6); // relative accuracy
	int dim, n;

	// Execute this for the four examples of dimensions listed in N:
	for (dim = 0; dim < 4; ++dim){
		n = N[dim];
		printf("Dimensions = n = %i \n", n);
		fprintf(f, "Dimensions = n = %i \n", n);
		printf("\n");
		fprintf(f, "\n");

		// Create pointers:
		double *a = calloc(n * n, sizeof(double));
		double *y = calloc(n, sizeof(double));
		double *x = calloc(n, sizeof(double));
		int *niter = calloc(1, sizeof(int)); 
		double *discreps = calloc(numit, sizeof(double)); 

		// Create the matrix A such that (A)i,i = 1/i^2 for i=1,2,...,n:
		int i, m;
		for (i = 1; i < n + 1; ++i){
			a[(i - 1) * n + (i - 1)] = 1 / pow(i, 2.0);
		}
		// Create vector y such that y(i) = 1 for i=1,2,...,n:
		for (m = 1; m < n + 1; ++m){
			y[m - 1] = 1.0;
		}
		
		/*
		// Check what I created:
		printf("A = \n");
		print_mat(a, n, n);
		printf("y = ");
		print_arr(y, n);
		printf("x = ");
		print_arr(x, n);
		*/

		// Call the subroutine dumb_solve:
		dumb_solve(a, y, n, eps, numit, x, niter, discreps);
		
		// Check results of dumb_solve:
		printf("Number of iterations actually performed = niter = %i\n", *niter);
		fprintf(f, "Number of iterations actually performed = niter = %i\n", *niter);
		printf("\n");
		fprintf(f, "\n");
		printf("Solution of the system after %i iterations: x = ", *niter); // n-dimensional array with coordinates of minimum after last iteration.
		fprintf(f, "Solution of the system after %i iterations: x = ", *niter);
		print_arr(f, x, n);
		printf("\n");
		fprintf(f, "\n");
		printf("Discrepancies after each iteration = "); // F(x) after each completed iteration
		fprintf(f, "Discrepancies after each iteration = ");
		print_arr(f, discreps, *niter); 
		printf("\n");
		fprintf(f, "\n");
		printf("_______________________________________________\n");
		fprintf(f, "_______________________________________________\n");
		
		// Free memory by destroying pointers:
		free(a); free(y); free(x); free(niter); free(discreps); 
	}
	fclose(f);
	return 0;
}