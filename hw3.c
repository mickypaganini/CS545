// --------------------------------------------------------------------//
// Michela Paganini                      							   //
// CS 545 - Introduction to Data Mining 							   // 
// Yale University                 								       //
// December 1, 2014                 							       //
// Homework 3                            							   //
// "NEAREST NEIGHBORS WITH ADAPTIVE QUAD-TREE"                         //
// GOAL: Find k nearest neighbors of each point in the plane		   //
// --------------------------------------------------------------------//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// ---------------------------------------//
//                Structures              //
// ---------------------------------------//
typedef struct MyArray {
    double data;
    int orig_pos; 
} MyArray;
// ---------------------------------------//
typedef struct point {
	double x;
	double y;
} point;
// ---------------------------------------//
typedef struct square {
	point bottomL;
	point bottomR;
	point topL;
	point topR;
} square;
// ---------------------------------------//
typedef struct table {
	square Square;
	int first;
	int last;
	int parent;
	int child[4];
} table;
// ---------------------------------------//
typedef struct circle {
	double radius;
	point origin;
} circle;
// ---------------------------------------//
//                Prototypes              //
// ---------------------------------------//
void seek(double *a, int n, int k, int *iz);
void seek_naive(double *a, int n, int k, int *iz);
void print_arr(FILE *f, int * arr, int n);
void print_mat(FILE *f, double * mat, int n_row, int n_col);
int compare(const void * a, const void * b);
void ordered_pairs(double *a, point *newarray, int n);
void give_birth(table *Table, point *A, int *index_array, int here, int next);
int isinthequadrant(point P, square S);
void swap(int *a, int *z);
int whosmydaddy(point P, table *T);
double radius(point P, square S);
double dist(point P, point C);
int overlap(circle K, square S);
// ---------------------------------------//
//                 Functions              //
// ---------------------------------------//
void seek(double *a, int n, int k, int *iz){
	// ################################################
	// Idea: 
	// TASK 1: create a table containing tree info.
	// 		   for each square, if it contains more thank k points, split it and store info in the table.
	//		   go thru all boxes until they all contain k or fewer elements.
	//		   fill the table with all the info for each box (location, number of elements, order of elements, parent, children)
	//
	// TASK 2: Loop thru all points; for a given point, find out in which square it is (aka go down th tree)
	//		   once we found the smallest box that contains it, find its parent box
	// 		   find the max distance of the point from any corner of the parent box and use that as a radius to draw a circle around the point
	//         consider all the boxes that intersect the circle: the points in those squares are the nearest neighbors candidates
	//         calculate distances of each point in those boxes to the point under investigation
	// ################################################

	// Error check:
	if (k > n - 1){
		printf("k needs to be greater than (n - 1).\n");
		return;
	}

	// Create A and index_array:
	point *A = malloc(n *sizeof(point)); // create array to store ordered points (x,y) to make my life simpler
	ordered_pairs(a, A, n); // fills A with data from a

	int *index_array = malloc(n *sizeof(int)); // create array to store indices to keep track of sorting
	int i;

	for (i = 0; i < n; ++i){ // fill index_array
		index_array[i] = i;
	}

	// Create table and first entry:
	table *Table = malloc(n *sizeof(table)); // make a huge table
	Table[0].Square.topL.x = 0; Table[0].Square.topL.y = 1;       // (0,1)
	Table[0].Square.topR.x = 1; Table[0].Square.topR.y = 1;       // (1,1)
	Table[0].Square.bottomL.x = 0; Table[0].Square.bottomL.y = 0; // (0,0)
	Table[0].Square.bottomR.x = 1; Table[0].Square.bottomR.y = 0; // (1,0)
	Table[0].first = 0; 
	Table[0].last = n - 1; 
	Table[0].parent = -999; // the first square is an orphan
	Table[0].child[0] = 0; Table[0].child[1] = 0; Table[0].child[2] = 0; Table[0].child[3] = 0;


	// ###########################################################################
	// CREATE TREE STRUCTURE (TASK 1):
	int here = 0; int next = here + 1;
	// Loop thru all the squares until we run out of squares
	while (here < next){

		// Split a square into 4 quadrants if it contains more than k elements:
		if ((Table[here].last - Table[here].first) > k){
			// printf("parto\n"); // make sure this executes
			give_birth(Table, A, index_array, here, next);
			next = next + 4; // if I split, next moves over by 4, if not it stays where it is
		}
		++here;

	}
	int last_box = next - 1;
	// ###########################################################################

	
	// ###########################################################################
	// GO DOWN PREVIOUSLY CREATED TREE STRUCTURE TO FIND k NEAREST NEIGHBORS (TASK 2): 
	// Loop thru all points:

	int j;
	for (j = 0; j < n; ++j){ // points loop
		
		// Find the parent box the point belongs to:
		int daddy_index = whosmydaddy(A[j], Table);
		
		// Find the distance of the point A[j] to the farthest corner of the parent and call that R:
		double R = radius(A[j], Table[daddy_index].Square);
		// Draw a circle of radius R centered at the point under consideration
		circle K;
		K.radius = R; K.origin = A[j]; 
		
		int box_number;
		int candidates_index[n]; // array containing the indices of possible nearest neighbors
		int idx = 0; // index to move thru array 'candidates_index'
		
		int yes_or_no[n]; // 1 if that point has already been identified as a candidate, 0 if it's not
		int yesno_i;
		for (yesno_i = 0; yesno_i < n; ++yesno_i){
			yes_or_no[yesno_i] = 0; // initialize it to an array of 0s
		}

		for (box_number = 0; box_number < (last_box + 1); ++box_number){ // box loop

			// Consider its parent: all points in the parent box are candidates
			if (box_number == daddy_index){

				int m;
				for (m = Table[box_number].first; m < (Table[box_number].last + 1); ++m){ // loop thru points in that box
					yes_or_no[index_array[m]] = 1; // identify this point as a candidate
					candidates_index[idx] = index_array[m];
					++idx;

				} 
			}
	
			// Consider all boxes that intersect the circle: all points in those boxes are candidates
			else if ((Table[box_number].child[0] == 0) && (Table[box_number].parent != daddy_index)){

				if (overlap(K, Table[box_number].Square)){

					int u;
					for (u = Table[box_number].first; u < (Table[box_number].last + 1); ++u){ // loop thru points in that box
						if (yes_or_no[index_array[u]] == 0){ // if it hasn't been identified yet
							candidates_index[idx] = index_array[u];
							++idx;
						}
					} 
				}
			}
			
		} // end loop thru boxes
		
		// Ok, so: now I have the 'candidates_index' array which contains all the indices of the candidate nearest neighbors
		// I can now go ahead an apply something similar to seek_naive to calculate the actual distances of these points and sort them.
		int y;
		MyArray *arr_of_distances = malloc(idx *sizeof(MyArray)); // array to store the distances of each point to point A[j]

		for (y = 0; y < idx; ++y){ // loop thru all candidates
			arr_of_distances[y].data = dist(A[j], A[candidates_index[y]]); // store their distances to A[j] 
			arr_of_distances[y].orig_pos = candidates_index[y]; 
		}
		
		// Reorder arr_of_distances:
		qsort(arr_of_distances, idx, sizeof(MyArray), compare);
		// Store the first k points into iz
		int z;
		for (z = 0; z < k; ++z){
			iz[j * k + z] = arr_of_distances[z + 1].orig_pos; // +1 because the first entry is 0 (= the distance from the point itself)
		}
		
		free(arr_of_distances); arr_of_distances = NULL;
		
	} // end loop thru points
	// ###########################################################################
	
	free(A); free(index_array); free(Table);
	A = NULL; index_array = NULL; Table = NULL;
}
// ---------------------------------------//
void seek_naive(double *a, int n, int k, int *iz){

	// Error check:
	if (k > n - 1){
		printf("k needs to be greater than (n - 1).\n");
		return;
		}

	int i, j, m, l, ind, index;
	double *distance = malloc(n * n *sizeof(double));
	MyArray *distancestopoint = malloc(n *sizeof(struct MyArray));

	// Create an nxn matrix of distances:
	for (i = 0; i < n; ++i){ // loop thru points to select one specific point

		for (j = i + 1; j < n; ++j){ // loop thru all other points
			distance[i * n + j] = sqrt((a[i * 2] - a[j * 2])*(a[i * 2] - a[j * 2]) + (a[i * 2 + 1] - a[j * 2 + 1])*(a[i * 2 + 1] - a[j * 2 + 1]));
			distance[j * n + i] = distance[i * n + j];
		}
	}

	// Throw every line of the matrix into an array called distancestopoint
	for (m = 0; m < n; ++m){ // loop thru rows aka points

		for (l = 0; l < n; ++l){ // loop thru elements of the row
			distancestopoint[l].data = distance[m * n + l]; // put distances into an array
			distancestopoint[l].orig_pos = l; // save the original order of indices
		}

		qsort(distancestopoint, n, sizeof(MyArray), compare); // sort in increasing order

		for (ind = 0; ind < k; ++ind){ // loop thru elements of newly ordered array
			iz[m * k + ind] = distancestopoint[ind + 1].orig_pos; // save the first k indices, excluding the first one because it's the point itself
		}
	}
	free(distance); free(distancestopoint);
	distance = NULL; distancestopoint = NULL;
}
// ---------------------------------------//
void print_arr(FILE *f, int * arr, int n){ 
// function that prints an array of size n 
	int i;
	printf("[ ");
	fprintf(f, "[ ");

	for (i = 0; i < (n - 1); ++i){
		printf("%i , ", arr[i]);
		fprintf(f, "%i , ", arr[i]);
	}
	printf("%i ]\n", arr[i]);
	fprintf(f, "%i ]\n", arr[i]);
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
int compare(const void * a, const void * b){
	// needed for qsort

	if ((*(MyArray*)a).data > (*(MyArray*)b).data)
   		return 1;

    else if ((*(MyArray*)a).data < (*(MyArray*)b).data)
    	return -1;

    else return 0;
}
// ---------------------------------------//
void ordered_pairs(double *a, point *newarray, int n){
	// it modifies the array 'newarray' by assigning it an x and y value from the array 'a'
	int i;
	for (i = 0; i < n; ++i){ // loop thru arrays 'a' and 'newarray'
		newarray[i].x = a[2 * i];
		newarray[i].y = a[2 * i + 1]; // y always follows x 
	}
}
// ---------------------------------------//
void give_birth(table *Table, point *A, int *index_array, int here, int next){
	// Generate children:
	int i;
	
	for (i = 0; i < 4; ++i){
		Table[here].child[i] = next + i; // fill the 'child' lines for the present square with location of children squares
		Table[next + i].parent = here; // fill children's 'parent' line with location of present square
	}
	
	// Assign coordinates to children quadrants:
	Table[next + 0].Square.topL.x = Table[here].Square.topL.x;                                              Table[next + 0].Square.topL.y = Table[here].Square.topL.y;
	Table[next + 0].Square.topR.x = (Table[here].Square.topR.x + Table[here].Square.topL.x)/ 2;             Table[next + 0].Square.topR.y = Table[here].Square.topR.y;
	Table[next + 0].Square.bottomL.x = Table[here].Square.bottomL.x; 	                                    Table[next + 0].Square.bottomL.y = (Table[here].Square.topL.y + Table[here].Square.bottomL.y) / 2;
	Table[next + 0].Square.bottomR.x = (Table[here].Square.bottomR.x + Table[here].Square.bottomL.x) / 2;   Table[next + 0].Square.bottomR.y = (Table[here].Square.bottomR.y + Table[here].Square.topR.y) / 2;

	Table[next + 1].Square.topL.x = (Table[here].Square.topL.x + Table[here].Square.topR.x) / 2;            Table[next + 1].Square.topL.y = Table[here].Square.topL.y;
	Table[next + 1].Square.topR.x = Table[here].Square.topR.x;                                              Table[next + 1].Square.topR.y = Table[here].Square.topR.y;
	Table[next + 1].Square.bottomL.x = (Table[here].Square.bottomL.x + Table[here].Square.bottomR.x) / 2;   Table[next + 1].Square.bottomL.y = (Table[here].Square.bottomL.y + Table[here].Square.topL.y) / 2;
	Table[next + 1].Square.bottomR.x = Table[here].Square.bottomR.x;                                        Table[next + 1].Square.bottomR.y = (Table[here].Square.bottomR.y + Table[here].Square.topR.y) / 2;

	Table[next + 2].Square.topL.x = Table[next + 0].Square.bottomL.x;                     					Table[next + 2].Square.topL.y = Table[next + 0].Square.bottomL.y;
	Table[next + 2].Square.topR.x = Table[next + 0].Square.bottomR.x; 										Table[next + 2].Square.topR.y = Table[next + 0].Square.bottomR.y;
	Table[next + 2].Square.bottomL.x = Table[here].Square.bottomL.x; 										Table[next + 2].Square.bottomL.y = Table[here].Square.bottomL.y; 
	Table[next + 2].Square.bottomR.x = (Table[here].Square.bottomR.x + Table[here].Square.bottomL.x) / 2;   Table[next + 2].Square.bottomR.y = Table[here].Square.bottomR.y;

	Table[next + 3].Square.topL.x = Table[next + 1].Square.bottomL.x; 										Table[next + 3].Square.topL.y = Table[next + 1].Square.bottomL.y;
	Table[next + 3].Square.topR.x = Table[next + 1].Square.bottomR.x; 										Table[next + 3].Square.topR.y = Table[next + 1].Square.bottomR.y;
	Table[next + 3].Square.bottomL.x = Table[next + 2].Square.bottomR.x;									Table[next + 3].Square.bottomL.y = Table[next + 2].Square.bottomR.y;
	Table[next + 3].Square.bottomR.x = Table[here].Square.bottomR.x; 										Table[next + 3].Square.bottomR.y = Table[here].Square.bottomR.y;
	
	// Assign points to each child quadrant:

	// Start with the 1st quadrant
	// Take the array of points and check if each point it in the 1st quadrant
	// If it is, increase the counter (aka move over one element) and leave point in its position in the 'index_array'
	// If it's not, swap it with the last non-checked element of 'index_array' and move over one elements from the tail --> come in from both sides
	// Do this until all elements of 'index_array' have been checked (aka when beginning>end). 
	// Then change quadrant and check the points that were not in the previous quadrants (aka from the new value of 'beginning')

	int quadrant;
	int first = Table[here].first; int last = Table[here].last; // initialize first and last, for shortcut

	for (quadrant = 1; quadrant < 5; ++quadrant){
		int beginning = first; int end = last;

		while (beginning <= end){

			if (isinthequadrant(A[index_array[beginning]], Table[next + quadrant - 1].Square)){
				beginning++;
			}

			else{
				swap(&index_array[beginning], &index_array[end]);
				end--;
			}
		}
		// Save info about elements in each child into the table:
		Table[next + quadrant - 1].first = first; Table[next + quadrant - 1].last = beginning - 1; 
		// Start next iteration from the points that have not been assigned to a quadrant yet
		first = beginning; 
	}

}
// --------------------------------------------//
int isinthequadrant(point P, square S){
	// function that checks whether a point is in the quadrant or not
	
	if ((P.x < S.topR.x) && (P.x > S.topL.x) && (P.y < S.topR.y) && (P.y > S.bottomR.y))
		return 1;

	else
		return 0;
}
// --------------------------------------------//
void swap(int *a, int *z){
	// function that swaps two integer entries in an array
	int temp = *a;
	*a = *z;
	*z = temp;
}
// --------------------------------------------//
int whosmydaddy(point P, table *T){
	// function that scans the entire tree to find the square, and therefore the parent that the point belongs to
	// we look for the square it's in, but ultimately only care about the index of the daddy
	// Start from the first entry of the table and look at its children.
	// If they exist, keep going down the tree. If they don't, stop:
	int box_index = T[0].child[0];
	int correct_box = 0; 

	while (box_index != 0){ // it would only be zero for boxes that don't have any children
		int k;

		for (k = 0; k < 4; ++k){
			if (isinthequadrant(P, T[box_index + k].Square)){
				correct_box = box_index + k;
				break;
			}
		}
		// Go one step down the tree:
		box_index = T[correct_box].child[0];
		//printf("box index = %i\n", box_index);
	} 
	return T[correct_box].parent;
}
// --------------------------------------------//
double radius(point P, square S){
	// function that finds distance of P to farthest corner of S
	double *dist_array = malloc(4 *sizeof(double));
	dist_array[0] = dist(P, S.topL);
	dist_array[1] = dist(P, S.topR);
	dist_array[2] = dist(P, S.bottomR);
	dist_array[3] = dist(P, S.bottomL);
	double maximum = dist_array[0];
	int l;

	for (l = 1; l < 4; l++){
    	if (dist_array[l] > maximum)
    		maximum  = dist_array[l];
    }
    free(dist_array); dist_array = NULL;
	return maximum;
}
// --------------------------------------------//
double dist(point P, point C){
	// function that calculates the distance between two points
	return sqrt((P.x - C.x) * (P.x - C.x) + (P.y - C.y) * (P.y - C.y));
}
// --------------------------------------------//
int overlap(circle K, square S){
	// function that checks whether the circle overlaps with the square
	if ((dist(S.topR, K.origin) <= K.radius) || (dist(S.topL, K.origin) <= K.radius) || (dist(S.bottomR, K.origin) <= K.radius) || (dist(S.bottomL, K.origin) <= K.radius))
		return 1;
	return 0;
}
// --------------------------------------------//

// --------------------------------------------//
//                  Main                       //
// 				  (Check)					   //
// --------------------------------------------// 
int main(){

	FILE *f;
	f = fopen("hw3.txt", "w"); // output file

	// Create and initialize quantities:
	int n = 1000000; // Number of points in the plane. It should be 1000000.
	int k = 100; // Number of nearest neighbors we are looking for. Pick whatever you want.
	double *a = malloc(n * 2 *sizeof(double)); // nx2 matrix
	int *iz = malloc(k * n *sizeof(int)); // array of size kn = array of indeces of the k NNs

	// Start by generating n points in a 2D plane by assigning random coordinates in x,y C [0,1]:
	int i, j;
	srand((unsigned)time(NULL)); // seed

	for (i = 0; i < n; ++i){

		for (j = 0; j < 2; ++j){
			a[i * 2 + j] = ((double)rand() / (double)RAND_MAX); // assign random coordinates to the points
		}
	}
	//printf("Matrix A: \n");
	//fprintf(f, "Matrix A: \n");
	//print_mat(f, a, n, 2); // check what I created

	// Use the naïve method to find the k nearest neighbors:
	seek_naive(a, n, k, iz);
	printf("From seek_naive:\n");
	fprintf(f, "From seek_naive:\n");
	printf("iz = ");
	fprintf(f, "iz = ");
	print_arr(f, iz, k * n);
	printf("\n");
	fprintf(f, "\n");
	
	// Reset iz:
	int ii;
	for (ii = 0; ii < (k * n); ++ii){
		iz[ii] = 0;
	}

	// Use the clever method to find the k nearest neighbors:
	seek(a, n, k, iz);
	printf("From seek:\n");
	fprintf(f, "From seek:\n");
	printf("iz = ");
	fprintf(f, "iz = ");
	print_arr(f, iz, k*n);

	free(a); free(iz);
	a = NULL; iz = NULL;
	fclose(f);
	return 0;
}