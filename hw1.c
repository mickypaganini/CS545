// ---------------------------------------//
// Michela Paganini                       //
// CS 545 - Introduction to Data Mining   // 
// Yale University                        //
// September 29, 2014                     //
// Homework 1                             //
// "SVD USING JACOBI ROTATION"            //
// ---------------------------------------//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// ---------------------------------------//
//                Prototypes              //
// ---------------------------------------//
void jacobi(double * a, int n, double * s, double * u, double * v);
void print_array(FILE *f, double * arr, int n);
void print_matrix(FILE *f, double * mat, int n_row, int n_col);
double* identity(int n);
double ATA_element(double *mat, int n, int i, int j);
double* jacobi_rotation(double *a, double *v, int i, int j, int n);
int sign(double x);
double* matrix_multiply(FILE *f, double *a, double *b, int rows_a, int cols_a, int rows_b, int cols_b);
double* transpose(double *mat, int n_row, int n_col);
double* sing_value(double *mat, int n);
int* new_index_order(double *arr, int n);
double* reorder(double *arr, int *index, int n);
double* order_columns(double *mat, int *index, int n);
double* findU(double *a, double *s, int n);
double* buildS(double *s, int n); 
void check_ortho(FILE *f, double *mat, int n, double eps);
// ---------------------------------------//
//                 Functions              //
// ---------------------------------------//
void jacobi(double * a, int n, double * s, double * u, double * v) {
	// ----------------------------------------------------------------------------//
	// DEFINITIONS:
	// Inputs:
	// a = nxn matrix to be diagonalized
	// n = dimensions
	// Outputs:
	// s = array of signular values in decreasing order = diag entries of U
	// u = left signular vectors = columns of U
	// v = right singular vectors = columns of V
	// a=usv^T
	// ----------------------------------------------------------------------------//
	FILE *f;
	f = fopen("hw1.txt", "w"); // output file
	printf("The original matrix A is: \n"); // write to screen
	fprintf(f, "The original matrix A is: \n"); // write to file
	print_matrix(f, a, n, n);

	// start by setting u and v equal to the n*n identity
	u = identity(n); v = identity(n);
	// pick most suitable epsilon to determine if off diag entry is already small 
	double eps = pow(10.0, -15), atemp, vtemp; 
	int exit_status = 1, i, j, row, itercount = 0;
	double *R = malloc(n * n *sizeof(double));

	while ((exit_status != 0) && (itercount < 10000)){
		exit_status = 0; // this breaks the loop if the value isn't changed
		//Cyclic-by-row-Jacobi: Sweep through the off diag elements of A row-wise
		for (i = 0; i < n - 1; ++i){ // rows
			for (j = i + 1; j < n; ++j){ // columns
					// only execute next line if A^T*A is not already diagonal
					if (fabs(ATA_element(a, n, i, j)) > (eps * sqrt(ATA_element(a, n, i, i) * ATA_element(a, n, j, j)))) { 
						// perform a Jacobi Rotation
						R = jacobi_rotation(a, v, i, j, n);
						// use R to modify a and v
						for (row = 0; row < n; ++row){ 
							atemp = a[row * n + i];
							a[row * n + i] = a[row * n + j] * R[j * n + i] + atemp * R[i * n + i];
							a[row * n + j] = a[row * n + j] * R[j * n + j] + atemp * R[i * n + j];
							vtemp = v[row * n + i];
							v[row * n + i] = v[row * n + j] * R[j * n + i] + vtemp * R[i * n + i];
							v[row * n + j] = v[row * n + j] * R[j * n + j] + vtemp * R[i * n + j];
						}
						exit_status += 1;
					}
			}
		}
		++itercount; // counts the iterations of the while loop
	}

	s = sing_value(a, n); // not in order
	int *indeces = new_index_order(s, n); // returns array of ordered indeces AND it reorders s!
	// --> now s is the array of singular values IN DECREASING ORDER
	a = order_columns(a, indeces, n); // reorder columns of a by decreasing singular value
	v = order_columns(v, indeces, n); // reorder columns of v by decreasing singular value
	double *vT = transpose(v, n, n); // transpose of v
	u = findU(a,s,n); // already in order of decreasing singular value
	double *S = buildS(s,n); // put singular values into diag matrix

	// outputs
	printf("The matrix A can be decomposed via SVD into A = US(V)^T, where: \n");
	fprintf(f, "The matrix A can be decomposed via SVD into A = US(V)^T, where: \n");
	printf("U is an orthogonal %i x %i matrix, containing the left singular vectors of A, given by: \n", n, n);
	fprintf(f, "U is an orthogonal %i x %i matrix, containing the left singular vectors of A, given by: \n", n, n);
	print_matrix(f,u,n,n);
	printf("S is a %i x %i diagonal matrix whose diagonal entries correspond to the singular values of A and are: \n", n, n);
	fprintf(f, "S is a %i x %i diagonal matrix whose diagonal entries correspond to the singular values of A and are: \n", n, n);
	print_array(f,s,n);
	printf("\n");
	fprintf(f, "\n");
	printf("(V)^T is an orthogonal %i x %i matrix given by: \n", n, n);
	fprintf(f, "(V)^T is an orthogonal %i x %i matrix given by: \n", n, n);
	print_matrix(f,vT,n,n);
	printf("Therefore, V is an orthogonal %i x %i matrix, containing the right singular vectors of A, given by: \n", n, n);
	fprintf(f, "Therefore, V is an orthogonal %i x %i matrix, containing the right singular vectors of A, given by: \n", n, n);
	print_matrix(f,v,n,n);

	// check that a=uS(v)^T
	printf("As a check, US(V)^T = A: \n");
	fprintf(f, "As a check, US(V)^T = A: \n");
	double *check = matrix_multiply(f,u,S,n,n,n,n); //uS
	check = matrix_multiply(f,check,vT,n,n,n,n); //uS(v)^T
	print_matrix(f,check,n,n);

	// check orthonormal
	printf("Checking if U is orthonormal: \n"); 
	fprintf(f, "Checking if U is orthonormal: \n");
	check_ortho(f,u,n,eps);
	printf("Checking if V is orthonormal: \n");
	fprintf(f, "Checking if V is orthonormal: \n");  
	check_ortho(f,v,n,eps);

	// free memory by destroying pointers
	free(R); free(indeces); free(S); free(check); free(vT); 
	fclose(f);
}
// ---------------------------------------//
void print_array(FILE *f, double * arr, int n){
	int i;
	for (i = 0; i < n; ++i){
		printf("%g ", arr[i]);
		fprintf(f, "%g ", arr[i]);
	}
	printf("\n");
	fprintf(f, "\n");
}
// ---------------------------------------//
void print_matrix(FILE *f, double * mat, int n_row, int n_col){
	int i, j;
	for (i = 0; i < n_row; ++i){
		for (j = 0; j < n_col; ++j){
			printf("%g ", mat[i * n_col + j]);
			fprintf(f, "%g ", mat[i * n_col + j]);
		}
		printf("\n");
		fprintf(f, "\n");
	}
	printf("\n");
	fprintf(f, "\n");
}
// ---------------------------------------//
double* identity(int n){
	double * I = malloc(n*n*sizeof(double));
	int k;
	for (k = 0; k < n; ++k){
		I[(n + 1) * k] = 1.0;
	}
	return I;
}
// ---------------------------------------//
double ATA_element(double *mat, int n, int i, int j){
	int row;
	double element_ij = 0;
	for(row = 0; row < n; ++row){ // loop thru rows of mat
		element_ij += mat[row * n + i] * mat[row * n + j];
	}
	return element_ij;
}
// ---------------------------------------//
double* jacobi_rotation(double *a, double *v, int i, int j, int n){
		double tau = (ATA_element(a, n, i, i) - ATA_element(a, n, j, j)) / (2 * ATA_element(a, n, i, j)); 
		double t = sign(tau) / (fabs(tau) + sqrt(1 + pow(tau, 2.0))); 
		double c = 1 / sqrt(1 + pow(t, 2.0)); // cosine
		double s = c * t; // sine=cosine*tangent
		double *R = identity(n);
		R[i*n+i] = c;
		R[j*n+j] = c;
		R[i*n+j] = -s;
		R[j*n+i] = s;
		return R;		
}
// ---------------------------------------//
int sign(double x){
	if (x<0)
		return -1;
	else if (x > 0)
		return +1;
	else
		return 0;
}
// ---------------------------------------//
double* matrix_multiply(FILE *f, double *a, double *b, int rows_a, int cols_a, int rows_b, int cols_b){
	if (cols_a != rows_b){
		printf("The number of columns of A must equal the number of rows of B.");
		fprintf(f, "The number of columns of A must equal the number of rows of B.");
		return NULL;
	}
	else{
		int row, col, index;
		double *product = malloc(rows_a * cols_b *sizeof(double));
		for (row = 0; row < rows_a; ++row){
			for (col = 0; col < cols_b; ++col){
				for (index = 0; index < cols_a; ++index){
					product[row * cols_b + col] += a[row * cols_a + index] * b[index * cols_b + col];
				}
			}
		}
		return product;
	}
}
// ---------------------------------------//
double* transpose(double *mat, int n_row, int n_col){
	double *trans = malloc(n_row * n_col *sizeof(double));
	int i,j;
	for (i = 0; i < n_col; ++i){
		for (j = 0; j < n_row; ++j){
			trans[i * n_row + j] = mat[j * n_col + i];
		}
	}
	return trans;
}
// ---------------------------------------//
double* sing_value(double *mat, int n){
	int i, j;
	double *s = malloc(n *sizeof(double));
	for (j = 0; j < n; ++j){ // loop thru columns
		double sum = 0;
		for (i = 0; i < n; ++i){ // loop thru rows
			sum += pow(mat[i * n + j],2); // sum of squares of elements of each column
		}
		s[j] = sqrt(sum); // sing value = sqrt of sum of squares of elements of each column
	}
	return s;
}
// ---------------------------------------//
int* new_index_order(double *arr, int n){
	int ind, k, l, temp_i;
	double temp;
	int* index = malloc(n *sizeof(int));
	// assign values to the entries of the array of indeces
	for (ind = 0; ind < n; ++ind){
		index[ind] = ind;
	}
	// swap until you get array containing the original indeces ordered from ind of max to ind of min entry
	for (k = 0; k < n; ++k){
		for (l = k + 1; l < n; ++l){
			if (arr[k] < arr[l]){
				temp = arr[k];
				arr[k] = arr[l];
				arr[l] = temp;
				temp_i = index[k];
				index[k] = index[l];
				index[l] = temp_i;
			}
		}
	}
	return index;
}
// ---------------------------------------//
double* reorder(double *arr, int *index, int n){
	int i,ind;
	double *ordered_arr = malloc(n*sizeof(double));
	for (i = 0; i < n; ++i){
		ind = index[i];
		ordered_arr[i] = arr[ind];
	}
	return ordered_arr;
}
// ---------------------------------------//
double* order_columns(double *mat, int *index, int n){
	double *new_mat = malloc(n * n*sizeof(double));
	int i, ind, k_row;
	for (i = 0; i < n; ++i){
		ind = index[i];
		for (k_row = 0; k_row < n; ++k_row){
			new_mat[k_row * n + i] = mat[k_row * n + ind];
		}
	}
	return new_mat;
}
// ---------------------------------------//
double* findU(double *a, double *s, int n){
	// {u_i} = {a_i}/(singvalue_i) where i labels columns
	double *u = malloc(n * n*sizeof(double));
	int j_col, i_row;
	for (j_col = 0; j_col < n; ++j_col){
		for (i_row = 0; i_row < n; ++i_row){
			u[i_row * n + j_col] = a[i_row * n + j_col] / s[j_col];
		}
	}
	return u;
}
// ---------------------------------------//
double* buildS(double *s, int n){
	double * S = malloc(n * n*sizeof(double));
	int k;
	for (k = 0; k < n; ++k){
		S[(n + 1) * k] = s[k]; // put sing values on diagonal
	}
	return S;
}
// ---------------------------------------//
void check_ortho(FILE *f, double *mat, int n, double eps){
	int checks = 0; // needs to be 2 for matrix to be orthonormal

	// 1) Check that all columns have length 1
	double *length_columns = sing_value(mat, n); //sqrt of sum of squared elements
	int i, non_norm_count = 0; 
	double compare;
	for (i = 0; i < n; ++i){
		compare = length_columns[i] - 1.0;
		if (fabs(compare) > sqrt(eps)){ // better comparison than epsilon itself
			printf("Column %i is not normalized. ", i);
			fprintf(f, "Column %i is not normalized. ", i);
			++non_norm_count;
		}
	}
	if (non_norm_count == 0){
		printf("The matrix is normalized.\n");
		fprintf(f, "The matrix is normalized.\n");
		++checks;
	}
	// 2) Check that all columns are orthogonal
	// take inner product of each column with all others
	int l, j, k, non_orth_count = 0;
	double inner_prod;
	for (l = 0; l < n; ++l){ // loop thru columns
		for (j = l + 1; j < n; ++j){ // loop thru other columns
			inner_prod = 0;
			for (k = 0; k < n; ++k){ // loop thru rows
				inner_prod += mat[l * n + k] * mat[j * n + k];
			}
			if (fabs(inner_prod) > sqrt(eps)){ 
				printf("Columns %i and %i are not orthogonal. \n", l, j);
				fprintf(f, "Columns %i and %i are not orthogonal. \n", l, j);
				++non_orth_count;
			}
		}
	}
	if (non_orth_count == 0){
		printf("All columns are orthogonal to one another. \n");
		fprintf(f, "All columns are orthogonal to one another. \n");
		++checks;
	}
	if (checks == 2){ // both 1) and 2) were true
		printf("--> The matrix is orthonormal. \n");
		fprintf(f, "--> The matrix is orthonormal. \n");
	}
	else {
		printf("--> The matrix is not orthonormal. \n");
		fprintf(f, "--> The matrix is not orthonormal. \n");
	}
	printf("\n");
	fprintf(f, "\n");
	free(length_columns);
}
// --------------------------------------------//
//                  Main                       //
// Only for checking purposes, can be deleted. //
// --------------------------------------------// 
int main(){
	// declare dimensions of matrix A
	int n_row = 6;
	int n_col = 6;
	// allocate mem for matrices I need
	double * a = malloc(n_row*n_col*sizeof(double));
	double * u = malloc(n_row*n_col*sizeof(double));
	double * v = malloc(n_row*n_col*sizeof(double));
	double *s = malloc(n_row*sizeof(double));
	// this just fills A with integers from 1 to 36
	int dummyindex;
	for (dummyindex = 0; dummyindex < n_row * n_col; ++dummyindex){
		a[dummyindex] = dummyindex + 1;
	}
	// call jacobi subroutine
	jacobi (a, n_col, s, u, v);
	free(a); free(u); free(v); free(s); 
	return 0;
}