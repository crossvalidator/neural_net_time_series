//#include "ancfis.h"
/*
This file contains code to calculate the checking (testing) error.
The checking error is calculated once per epoch.

There are 2 functions:
1. epoch_checking_error
	2. checking_error_measure
*/

/***********************************************
Description:
This function is responsible for calculating the checking data output values.
It then calls a function to calculate various error measures and generate an output file.

Inputs:
**data_matrix - contains the checking data from the checking data file
anfis_output - the predicted output for all checking data pairs
chk_data_n - the number of checking data pairs
error_index[] - stores the error results
trn_data_n - the number of training data pairs
epoch - the current epoch

Outputs:
error_index[] is considered to be the output.

Called from:
mexFunction in ancfismex.c
************************************************/
void epoch_checking_error(double **data_matrix, int chk_data_n, double *error_index_n, double *error_index_un, int trn_data_n, double **ancfis_output, int epoch, double minmaxc[]) {
	
	double total_squared_error_n=0;
	double total_squared_error_un=0;
	double RMSE=0;
	double diff0=0,diff1=0, diff2=0;
	int i,j;
	FILE *fp;
	FILE * open_file();

	void put_input_data();
	void calculate_output();
	void checking_error_measure();

	
	fp = open_file("tracking.txt", "a");
	/*if((fp = fopen("tracking.txt", "w")) == NULL)
	{
		printf("Cannot open 'data_file'.\n");
	}*/
	
	fprintf(fp, "\n=====================\nepoch: %d\n", epoch+1);

	if(chk_data_n == 0) {
		printf("No data in given data matrix!\n");
	}

	for(j = 0; j < chk_data_n; j++) {
		put_input_data(node_p,j, data_matrix);
		calculate_output(In_n, Node_n, j);
		for (i=0; i<Out_n; i++)
		{
			ancfis_output[j][i] = node_p[Node_n - (Out_n-i)]->value->real;
		}
	}
/**********************/
	/*for(i = 0; i < chk_data_n; i++)
	{
		for(j = 0; j < Out_n; j++)
		{	
			fprintf(fp, "%d \t	%lf \t   %lf\n",j, data_matrix[i][(j+1)*In_vect_n+j], ancfis_output[i][j]);
			diff0 = (data_matrix[i][(j+1)*In_vect_n+j] - ancfis_output[i][j]);
			diff0 = diff0*diff0;
			diff1 += diff0;
		}
		total_squared_error += (diff1/Out_n);
		error_index[j]=sqrt(total_squared_error);
	}*/

	for(i = 0; i < Out_n; i++)
	{
		for(j = 0; j < chk_data_n; j++)
		{	
			fprintf(fp, "%d \t	%lf \t   %lf\n",j, data_matrix[j][(i+1)*In_vect_n+i], ancfis_output[j][i]);
			diff0 = (data_matrix[j][((i+1)*(In_vect_n))+i] - ancfis_output[j][i]);
			diff0 = diff0*diff0;
			diff1 += diff0;
		}
		diff2 = diff1 * ((minmaxc[(2*i) +1] - minmaxc[2*i])*(minmaxc[(2*i) +1] - minmaxc[2*i]));
		total_squared_error_n += diff1;
		total_squared_error_un += diff2;
		diff1= (double)(diff1/chk_data_n);
		diff2= (double)(diff2/chk_data_n);
		error_index_n[i]=sqrt(diff1);
		error_index_un[i]=sqrt(diff2);
		diff1=0;
		
	}
	//printf("testing 5 \n");	
	//printf("testing 6 \n");
	fclose(fp);
	
	RMSE = sqrt(total_squared_error_n/(Out_n*chk_data_n));
	error_index_n[Out_n]=RMSE;
	RMSE = sqrt(total_squared_error_un/(Out_n*chk_data_n));
	error_index_un[Out_n]=RMSE;
	//checking_error_measure(data_matrix, ancfis_output, chk_data_n, error_index, epoch);
//	checking_error_measure(data_matrix, ancfis_output, chk_data_n,Out_n, error_index);
}


