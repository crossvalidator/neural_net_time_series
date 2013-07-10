//#include "ancfis.h"
/*
This file contains various 'library' function.  This includes
memory allocation/deallocation functions.

There are 9 functions:
1. create_matrix
2. create_array
3. open_file
4. free_matrix
5. free_array
6. print_matrix
7. print_array
8. write_matrix
9. write_array
*/

/***********************************************
Description:
Allocates memory for a matrix (2d array) using
Matlab's memory allocation functions.

Inputs:
row_n - the number of rows in the matrix
col_n - the number of columns in the matrix
element_size - the size of each element of the matrix

Outputs:
Allocates memory and returns the matrix

Called from:
mexFunction in ancfismex.c
initpara in initpara.c
kalman in kalman.c
************************************************/
char **create_matrix(int row_n, int col_n, size_t element_size) {
//	char **create_matrix(int row_n, int col_n, int element_size) {

	char **matrix;
	int i;

	matrix = calloc(row_n, sizeof(char *));
	if(matrix == NULL)
		printf("Error in create_matrix1!\n");
	for(i = 0; i < row_n; i++) {
		matrix[i] = calloc(col_n, element_size);
		if(matrix[i] == NULL)
			printf("Error in create_matrix2!\n");
	}

	return(matrix);
}

/***********************************************
Description:
Create an array with length 'array_size' and element size 'element_size'
using Matlab's memory allocation functions

Inputs:
array_size - the desired length of the array
element_size - the size of each element

Outputs:
Allocates the memory for the array and returns it

Called from:
mexFunction in ancfismex.c
kalman in kalman.c
************************************************/
char *create_array(int array_size, int element_size) {
	char *array;

	array = calloc(array_size, element_size);
	if(array == NULL)
		printf("Error in create_array!\n");

	return(array);
}

/***********************************************
Description:
A friendly interface to fopen()

Inputs:
*file - the file name
*mode - read, write, etc.

Outputs:
A file pointer to the opened file.

Called from:
mexFunction in ancfismex.c
write_all in debug.c
write_array and write_matrix in lib.c
write_parameter in output.c
************************************************/
FILE *open_file(char *file, char *mode) {
	FILE *fp, *fopen();

	if((fp = fopen(file, mode)) == NULL) {
		printf("Cannot open file");
	}

	return(fp);
}

/***********************************************
Description:
Frees a matrix created by create_matrix

Inputs:
**matrix - the matrix to be freed
row_n - the number of rows in the matrix

Outputs:
N/A

Called from:
N/A
************************************************/
void free_matrix(char **matrix, int row_n) {
	int i;

	for (i = 0; i < row_n; i++)
		if (matrix[i] == NULL) {
			//do nothing
		} else {
			free(matrix[i]);
			matrix[i] = NULL;
		}

	if (matrix == NULL) {
		//do nothing
	} else {
		free(matrix);
		matrix = NULL;
	}
}

/***********************************************
Description:
Frees a matrix created by create_array

Inputs:
array - the array to be freed

Outputs:
N/A

Called from:
N/A
************************************************/
void free_array(char *array) {
	if (array == NULL) {
		//do nothing
	} else {
		free(array);
		array = NULL;
	}
}

/***********************************************
Description:
Prints out a matrix to the screen.

Inputs:
matrix - the matrix to be printed
row_n - the number of rows in the matrix
col_n the number of columns in the matrix

Outputs:
N/A

Called from:
N/A
************************************************/
void print_matrix(int **matrix, int row_n, int col_n) {
	int i, j;
	for (i = 0; i < row_n; i++) {
		for (j = 0; j < col_n; j++)
			printf("%d ", matrix[i][j]);
		printf("\n");
	}
}

/***********************************************
Description:
Print out an array to the screen.

Inputs:
*array - the array to be printed
size - the number of elements in the array

Outputs:
N/A

Called from:
N/A
************************************************/
void print_array(double *array, int size) {
	int i;
	for (i = 0; i < size; i++)
		printf("%lf\n", array[i]);
}

/***********************************************
Description:
Writes a matrix to a file.

Inputs:
**matrix - the matrix to be written to a file
row_n - the number of rows in the matrix
col_n - the number of columns in the matrix
file_name - the name of the file to be written

Outputs:
Creates a file.

Called from:
N/A
************************************************/
void write_matrix(double **matrix, int row_n, int col_n, char *file_name) {
	FILE *fp;
	int i, j;

	fp = open_file(file_name, "w");
	for (i = 0; i < row_n; i++) {
		for (j = 0; j < col_n; j++)
			fprintf(fp, "%.11f ", matrix[i][j]);
		fprintf(fp, "\n");
	}
}

/***********************************************
Description:
Write an array to a file.

Inputs:
*array - the array to be written to a file
size - number of elements in the array
*file_name - the name of the file to be written

Outputs:
Creates a file.

Called from:
mexFunction in ancfismex.c
************************************************/
void write_array(double *array, int size, char *file_name) {
	FILE *fp;
	int i;

	fp = open_file(file_name, "w");
	for (i = 0; i < size; i++)
		fprintf(fp, "%d %.11f\n", i+1, array[i]);
	fclose(fp);
}

