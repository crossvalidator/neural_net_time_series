//#include "ancfis.h"
//#include "mex.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
/*
This program takes in a time series data file with each entry in its own row
(the input data file should have 1 column) and formats it for use in ANCFIS.

eg for a file with (1,2,3,4,5,6,7,8,9,10)

ip_vector_length = 4 with all else at 1 will give:
1 2 3 4|5
2 3 4 5|6
3 4 5 6|7
4 5 6 7|8
5 6 7 8|9
6 7 8 9|10

ip_vector_length = 4 and
ip_spacing = 2 with all else at 1 will give:
1 3 5 7|8
2 4 6 8|9
3 5 7 9|10

ip_vector_length = 4 and
ip_op_space = 2 with all else at 1 will give:
1 2 3 4|6
2 3 4 5|7
3 4 5 6|8
4 5 6 7|9
5 6 7 8|10

ip_vector_length = 2 and
ip_spacing = 2 and
ip_op_space = 2 and
op_vector_length = 2 with all else at 1 will give:
1 3|5 6
2 4|6 7
3 5|7 8
4 6|8 9
5 7|9 10

ip_vector_length = 3 and
op_vector_length = 3 and
op_spacing = 2 with all else at 1 will give:
1 2 3|4 6 8
2 3 4|5 7 9
3 4 5|6 8 10

ip_vector_length = 3 and
pair_spacing = 2 with all else at 1 will give
1 2 3|4
3 4 5|6
5 6 7|8
7 8 9|10

ip_vector_length = 2 and
ip_spacing = 2 and
ip_op_space = 3 and
op_vector length = 2 and
op_spacing = 2 and
pair_spacing = 3 gives
1 3|6 8
*/
/***********************************************
Description:


Inputs:


Outputs:


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




//void main(void)
int main()
{
	int my_array[10];
	int i, j, k, total_span, num_iterations_training, num_iterations_checking;
	
	int number_of_rows_training, number_of_rows_checking, ip_vector_length, ip_spacing, ip_op_spacing, op_vector_length, op_spacing, pair_spacing, normalize, In_n;
	FILE *input_file, *output_file, *fpp;
	double **trn_data, **chk_data, tmp;// mint, maxt, minc, maxc;

/****************************************************************
	//double dum_trn_data
	//string prhs[]={op_file_name, number_of_rows, ip_vector_length, ip_spacing, ip_op_spacing, op_vector_length, op_spacing, pair_spacing, ip_normalize, op_normalize};
   // if (nrhs != 11) {
	//     printf("Usage: file_name, op_file_name, number_of_rows_in_file, input_vector_length, input_spacing, space_between_last_input_and_first_output, output_vector_length, output_spacing, space_between_data_pairs ip_normalize op_normalize\n");
	//}

*****************************************************************/


	//int normalization;

	char op_file_name_training []= "data.trn";
	char op_file_name_checking []= "data.chk";
	char file_name_training []= "motel.dat.trn";
	char file_name_checking []= "motel.dat.chk";


	//file_name="1909-1979.txt";
	//op_file_name="data.chk";

	//file_name = (char[])calloc(17, sizeof(char)); &&&&&&&&&& i have chnaged this step to the folwoing
           //input_file= (FILE *)calloc(17, sizeof(char));

/*	status = gets(prhs->In_n, file_name, buflen);
	if (status != 0) {
    	printf("Not enough space. String is truncated.");
	}
*/
	
	//op_file_name = (char[] *)calloc(17, sizeof(char));  &&&&&&&&&& i have chnaged this step to the folwoing
      // output_file= (FILE *)calloc(17, sizeof(char));

/*****************************************************************************************
	status = mxGetString(prhs[1], op_file_name, buflen);
	if (status != 0) {
    	printf("Not enough space. String is truncated.");
	}
***************************************************************************************************/


// reading the vital parametrs from the datapara.txt file
	if((fpp = fopen("datapara.txt", "r")) == NULL)
	{
		printf("Cannot open 'parameter_file'.\n");
	}

	for(j = 0; j < 10; j++)
	{
		if(fscanf(fpp, "%lf", &tmp) != EOF) 
		{
				my_array[j] = tmp;
				
		} else 		{
			printf("Not enough data in 'input_parameter'!");
		}
	}

	fclose(fpp);
 
// assigning values to all the vital parametrs

   	number_of_rows_training = my_array[0];
	number_of_rows_checking = my_array[1];
	ip_vector_length = my_array[2];
	ip_spacing = my_array[3];
	ip_op_spacing = my_array[4];
	op_vector_length= my_array[5];
	op_spacing= my_array[6];
	pair_spacing= my_array[7];
	normalize= my_array[8];
	In_n = my_array[9];
/**************************************************************************************
	//TRAINING
	number_of_rows =221; //min 1
	ip_vector_length =12 ; //min 1
	//TESTING
	//number_of_rows =71; //min 1
	//ip_vector_length =12 ; //min 1
	
	ip_spacing =1; //min 1
	ip_op_spacing =1; //min 1
	op_vector_length =1; //min 1
	op_spacing =1; //note this is irrelevant if op_vector_length is 1
	pair_spacing =1; //min 1
	ip_normalize =1; //0 or 1
	op_normalize =1; //0 or 1
	//min = DBL_MAX;
	//max = DBL_MIN;
*********************************************************************************/

 double mint[In_n], maxt[In_n], minc[In_n], maxc[In_n];

	if (number_of_rows_training < 1 || number_of_rows_checking < 1 || ip_vector_length < 1 || ip_spacing < 1 || ip_op_spacing < 1 || op_vector_length < 1 || op_spacing < 1 || pair_spacing < 1 || In_n < 1) 
	{
		printf("Spacing parameters cannot be < 1.\n");
	}


	//store the time series data in an array
	//trn_data = (double **)calloc(number_of_rows_training, In_n, sizeof(double));
	trn_data =(double **)create_matrix(number_of_rows_training, In_n, sizeof(double)); 
	//double dum_trn_data[number_of_rows];
	chk_data = (double **)create_matrix(number_of_rows_checking, In_n, sizeof(double)); 

// getting the training data from the training data file in trn_data 
	if ((input_file = fopen(file_name_training, "r")) == NULL) {
		//fprintf("Cannot open datafile.\n");
		printf("Cannot open datafile.\n");
	}

	//free(file_name_training);

	for (i = 0 ; i < number_of_rows_training ; i++) 
	{
		for(j=0 ; j < In_n ; j++)
			{
				fscanf(input_file, "%lf", &tmp);
				trn_data[i][j] = tmp;
			}
	}

	fclose(input_file);

// getting the data from cheking data file to chk_data
	if ((input_file = fopen(file_name_checking, "r")) == NULL) {
		//fprintf("Cannot open datafile.\n");
		printf("Cannot open datafile.\n");
	}

	//free(file_name_checking);

	
	for (i = 0 ; i < number_of_rows_checking ; i++) 
	{
		for(j=0 ; j < In_n ; j++)
			{
				fscanf(input_file, "%lf", &tmp);
				chk_data[i][j] = tmp;
			}
	}


	fclose(input_file);

	for(i=0; i<In_n; i++)
	{	
	mint[i]=trn_data[0][i];
	maxt[i]=trn_data[number_of_rows_training -1][i];

	minc[i]=chk_data[0][i];
	maxc[i]=chk_data[number_of_rows_checking -1][i];
	}


// finding the maximum and the minimum of the training data pairs and then normalizing
	for(j=0; j<In_n; j++)
	{	
	for (i = 0 ; i < number_of_rows_training ; i++) {
		if (trn_data[i][j] < mint[j]) {
			mint[j] = trn_data[i][j];
		}

		if (trn_data[i][j] > maxt[j]) {
			maxt[j] = trn_data[i][j];
		}
	}
	}

	for(j=0; j<In_n; j++)
	{
	if((mint[j] != maxt[j])&& (normalize ==1))
	{
	for (i = 0 ; i < number_of_rows_training ; i++) 
	{
		trn_data[i][j]=(trn_data[i][j]-mint[j])/(maxt[j]-mint[j]);
	}
	}
	}


// finding the minimum and maximum of the checking data pairs and then normalizing
	for(j=0; j<In_n; j++)
	{		
	for (i = 0 ; i < number_of_rows_checking ; i++) {
		if (chk_data[i][j] < minc[j]) {
			minc[j] = chk_data[i][j];
		}

		if (chk_data[i][j] > maxc[j]) {
			maxc[j] = chk_data[i][j];
		}
	}
	}


	if ((output_file = fopen("minmax.txt", "w")) == NULL) {
		printf("Cannot write training data.\n");
	}
	
	for(i=0; i<In_n; i++)
	{
	
	fprintf(output_file, "%lf \n",minc[i]);
	fprintf(output_file, "%lf \n", maxc[i]);
	}
	
	fclose(output_file);

	for(j=0; j<In_n; j++)
	{	
	if((minc[j] != maxc[j])&& (normalize ==1))
	{
	for (i = 0 ; i < number_of_rows_checking ; i++) 
	{
		//chk_data[i][j]=(chk_data[i][j]-minc[j]);
		//printf("checking 1st line %lf \t", chk_data[i][j]);
		chk_data[i][j]=(chk_data[i][j]-minc[j])/(maxc[j]-minc[j]);
		//printf("checking 1st line %lf \t", chk_data[i][j]);
	}
	}
		//printf("\n");
	}

	
/*********************************************************************************************
	//printf("%f, %f", min, max);
/*
	
	if (ip_normalize == 1) {
		if (min < 0) {
			for (i = 0 ; i < number_of_rows ; i++) {
				trn_data[i] += fabs(min);
			}

			max += fabs(min);
		}

		//there could be problems if max = min
		//ex: all the data is the same, 0 0 0 0 0 or -2 -2 -2 -2
		if(max != min) {
			for (i = 0 ; i < number_of_rows ; i++) {
			if(min<0)
			     trn_data[i] /= max;
            else
				trn_data[i] = (trn_data[i]-min)/(max-min);
				//printf("%lf\n",trn_data[i]);
			}
		} else {
			printf("ERR in normalize");
		}
	}


	*/

/*
	//calculate the number of iterations

	//need to get max span of 1 iteration
	//span of input is (ip_vector_length-1)*ip_spacing
	//span b/w input and output is ip_op_spacing
	//span b/w outputs is (op_vector_length-1)*op_spacing
	//total_span = those 3 spans + 1

	//num_iterations = 1 + (number_of_rows-total_span)/pair_spacing

********************************************************************************************************/
	total_span = 1 + (ip_vector_length - 1)*ip_spacing + ip_op_spacing + (op_vector_length - 1)*op_spacing;
	num_iterations_training = 1 + floor((number_of_rows_training - total_span)/pair_spacing);
	num_iterations_checking = 1 + floor((number_of_rows_checking - total_span)/pair_spacing);

	if ((output_file = fopen(op_file_name_training, "w")) == NULL) {
		printf("Cannot write training data.\n");
	}

	//free(op_file_name_training);
	
	for (i = 0 ; i < num_iterations_training ; i++) {
		for(k=0; k<In_n; k++)
		{
		//the input vector
		for (j = 0 ; j < ip_vector_length ; j++) {
			fprintf(output_file, "%4.10lf ", trn_data[i*pair_spacing+j*ip_spacing][k]);
		}
		//the output vector
		for(j = 0 ; j < op_vector_length ; j++) {
			/*
				if(op_normalize == 0 && ip_normalize == 1) { //undo the normalization
				tmp = trn_data[i*pair_spacing+(ip_vector_length - 1)*ip_spacing + ip_op_spacing + j*op_spacing]*(max-min)+min;
			} else {
				tmp = trn_data[i*pair_spacing+(ip_vector_length - 1)*ip_spacing + ip_op_spacing + j*op_spacing];
			}
			*/

			tmp = trn_data[i*pair_spacing+(ip_vector_length - 1)*ip_spacing + ip_op_spacing + j*op_spacing][k];

			if(j != op_vector_length - 1) {
				fprintf(output_file, "%4.10lf ", tmp);
			} else {
				fprintf(output_file, "%4.10lf", tmp);
			}
		}

		if((i != num_iterations_training - 1)||(k != In_n-1)) {
			fprintf(output_file, "\n");
								      }
		}
		}


	fclose(output_file);

	if ((output_file = fopen(op_file_name_checking, "w")) == NULL) {
		printf("Cannot write training data.\n");
	}

	//free(op_file_name_checking);

	for (i = 0 ; i < num_iterations_checking ; i++) {
		for(k=0; k<In_n; k++)
		{
		//the input vector
		for (j = 0 ; j < ip_vector_length ; j++) {
			fprintf(output_file, "%4.10lf ", chk_data[i*pair_spacing+j*ip_spacing][k]);
		}
		//the output vector
		for(j = 0 ; j < op_vector_length ; j++) {
			/*
				if(op_normalize == 0 && ip_normalize == 1) { //undo the normalization
				tmp = trn_data[i*pair_spacing+(ip_vector_length - 1)*ip_spacing + ip_op_spacing + j*op_spacing]*(max-min)+min;
			} else {
				tmp = trn_data[i*pair_spacing+(ip_vector_length - 1)*ip_spacing + ip_op_spacing + j*op_spacing];
			}
			*/

			tmp = chk_data[i*pair_spacing+(ip_vector_length - 1)*ip_spacing + ip_op_spacing + j*op_spacing][k];

			if(j != op_vector_length - 1) {
				fprintf(output_file, "%4.10lf ", tmp);
			} else {
				fprintf(output_file, "%4.10lf", tmp);
			}
		}

		if((i != num_iterations_checking - 1)|| (k!= In_n -1)) {
			fprintf(output_file, "\n");
		}
	}
	}
	fclose(output_file);

	

	


	
	free(trn_data);
	free(chk_data);
return(0);
	//plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	//matlab_output = mxGetPr(plhs[0]);
	//*matlab_output = num_iterations;
}

