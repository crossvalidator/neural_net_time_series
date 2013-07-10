//#include "ancfis.h"
#include <stdio.h>
#include <stdlib.h>
/*
This file contains functions that may be useful for debugging ANCFIS

There are 9 functions:
1. debug_anfis
2. print_fan_in
3. print_fan_out
4. print_parameter
5. print_all_parameter
6. pause
7. pause2
8. write_all
9. get_file_name
*/

/***********************************************
Description:
Pauses execution, prompts user to input a node index.
The fan-in and fan-out nodes, along with other information is printed
to screen for that node.

Inputs:
N/A

Outputs:

Called from:
Wherever the user would like to place it.
************************************************/
void debug_anfis() {
	void print_fan_in(), print_fan_out(), print_parameter();
	int index;
//	mxArray *input, *str;


	while(1) {
	//	mexCallMATLAB(1, &input, 1, &str, "input");
		//index = (int) mxGetScalar(input);
		index=1;

		if(index == -1) {
			break;
		}

		if(index < 0 || index > Node_n - 1) {
			printf("Invalid node index!\n");
			continue;
		}

		printf("\n");
		printf("current node: %d [value = (%f, %f), de_do = (%f, %f), layer = %d, f_index = %d]\n",
		index, node_p[index]->value->real,
		node_p[index]->value->imag,
		node_p[index]->de_do->real,
		node_p[index]->de_do->imag,
		node_p[index]->layer,
		node_p[index]->function_index);
		print_fan_in(index);
		print_fan_out(index);
		print_parameter(index);
		printf("\n");
	};

	//mxDestroyArray(input);
	//mxDestroyArray(str);
}

/***********************************************
Description:
Prints the fan in list of node_p[index].

Inputs:
index - the desired index

Outputs:
N/A

Called from:
Wherever the user would like to place it.
debug_anfis in debug.c
************************************************/
void print_fan_in(int index) {
	NODE_LIST_T *p;

	if(node_p[index]->fan_in == NULL) {
		printf("No fan-in nodes!\n");
		return;
	}

	printf("fan-in nodes: ");

	for(p = node_p[index]->fan_in; p != NULL; p = p->next)
		printf("%d ", p->content->index);

	printf("\n");
}

/***********************************************
Description:
Prints the fan out list of node_p[index].

Inputs:
index - the desired index

Outputs:
N/A

Called from:
Wherever the user would like to place it.
debug_anfis in debug.c
************************************************/
void print_fan_out(int index) {
	NODE_LIST_T *p;

	if(node_p[index]->fan_out == NULL) {
		printf("No fan-out nodes!\n");
		return;
	}

	printf("fan-out nodes: ");

	for(p = node_p[index]->fan_out; p != NULL; p = p->next)
		printf("%d ", p->content->index);

	printf("\n");
}

/***********************************************
Description:
Prints the parameter information at node_p[index].

Inputs:
index - the desired index

Outputs:
N/A

Called from:
Wherever the user would like to place it.
debug_anfis in debug.c
************************************************/
void print_parameter(int index) {
	PARAMETER_LIST_T *p;

	if(node_p[index]->parameter == NULL) {
		printf("No parameters for this node!\n");
		return;
	}

	printf("parameter values: (");

	for(p = node_p[index]->parameter; p != NULL; p = p->next)
		printf("%f ",p->content);

	printf(")\n");
}

/***********************************************
Description:
Print ALL parameter information as well as de_dp

Inputs:
N/A

Outputs:
N/A

Called from:
Wherever the user would like to place it.
************************************************/
void print_all_parameter() {
	PARAMETER_LIST_T *p;
	int i;

	for(i = 0; i < Node_n; i++) {
		if (node_p[i]->parameter == NULL)
			continue;

		for(p = node_p[i]->parameter; p != NULL; p = p->next) {
			printf("%4.5lf ", p->content);
			printf("(%4.5lf, %%4.5lf) ", p->de_dp->real, p->de_dp->imag);
		}

		printf("\n");
	}
}


/***********************************************
Description:
Simply pauses the execution of the code until the user inputs
from the keyboard.

Insert this statement into the code and recompile to use it.

Inputs:
N/A

Outputs:
N/A

Called from:
Wherever the user would like to place it.
************************************************/
/*void pause() {
	mxArray *input, *str;
	str = mxCreateString("Please any key to continue");
	mexCallMATLAB(1, &input, 1, &str, "input");
	mxDestroyArray(input);
	mxDestroyArray(str);
}
*/
/***********************************************
Description:
Pauses the execution of the code with a custom message until
the user inputs frmo the keyboard.

Inputs:
*msg - a string message that will be displayed on the screen
	when it is paused.

Outputs:
N/A

Called from:
Wherever the user would like to place it.
************************************************/
/*void pause2(char *msg) {
	mxArray *input, *str;
	str = mxCreateString(msg);
	mexCallMATLAB(1, &input, 1, &str, "input");
	mxDestroyArray(input);
	mxDestroyArray(str);
}*/

/***********************************************
Description:
APPENDS all parameters, nodes & their values to a file 'log.txt'.
This can be used to verify the contents of all of the nodes AFTER execution.

Inputs:
epoch_n - the current epoch
training_data_n - the number of training data pairs
*msg - a string in case the user would like to give a description of what
	is going to be outputted into the file.

Outputs:
An output file is created in the current working directory called 'log.txt'.

Called from:
Wherever the user would like to place it.
************************************************/
void write_all(int epoch_n, int training_data_n, char *msg) {
	FILE *fp;
	PARAMETER_LIST_T *p;
	int i,k,h;

	FILE * open_file();

	fp = open_file("log.txt", "a");
	fprintf(fp, msg);
	fprintf(fp, "Epoch: %d Pair#: %d\n", epoch_n, training_data_n);

	for(i = 0; i < Node_n; i++) {
		if(i<In_n){
			fprintf(fp,"node:%d   output=", i);
        		for(k = 0; k < In_n; k++) {
			double *trn_data = node_p[k]->input_vector;
			for(h = 0 ; h < In_vect_n ; h++) 
				fprintf(fp, "%4.5lf   ",trn_data[h]);}
        fprintf(fp,"de_do = (%4.5lf, %4.5lf), layer = %d, f_index = %d]\n",
		node_p[i]->de_do->real, node_p[i]->de_do->imag,
		node_p[i]->layer,
		node_p[i]->function_index);
		}
		else
		{
		fprintf(fp, "node: %d [output = (%4.5lf, %4.5lf), de_do = (%4.5lf, %4.5lf), layer = %d, f_index = %d]\n",
		i,
		node_p[i]->value->real, node_p[i]->value->imag,
		node_p[i]->de_do->real, node_p[i]->de_do->imag,
		node_p[i]->layer,
		node_p[i]->function_index);}
		if(node_p[i]->parameter == NULL) {
			continue;
		} else {
			fprintf(fp, "Parameters:\n");	
		}

		for(p = node_p[i]->parameter; p != NULL; p = p->next) {
			fprintf(fp, "%4.5lf ", p->content);
			fprintf(fp, "de_dp = (%4.5lf, %4.5lf)\n", p->de_dp->real, p->de_dp->imag);
		}
	}

	fprintf(fp, "=======\n");

	fclose(fp);
}

/***********************************************
Description:
Generates a string based on the:
stepsize, stepsize increase rate, stepsize decrease rate,
number of membership functions and input vector length.

This string is intended to be a file name for cases where
a user is trying many different combinations of parameters to find
the best training/checking error.

Inputs:
ss - stepsize
ss_d - stepsize decrease rate
ss_i - stepsize increase rate

Outputs:
The output is a string. Note that decimal only print out to 2 places.

Ex: ss = 0.2, ss_d = 0.9, ss_i = 1.15, mf_n = 3, vect_len = 42 gives
"0.20_0.90_1.15_42_3"

Called by:
Wherever the user would like to place it.
************************************************/
char *get_file_name(double ss, double ss_d, double ss_i) {
	char ssize[10];
	char ssize_i[10];
	char ssize_d[10];
	char vect_len[10];
	char mf_n[10];
	char *file_name;

	sprintf(ssize, "%.2f", ss);
	sprintf(ssize_d, "%.2f", ss_d);
	sprintf(ssize_i, "%.2f", ss_i);
	sprintf(vect_len, "%d", In_vect_n);
	sprintf(mf_n, "%d", Mf_n);

	file_name = calloc(strlen(ssize) + strlen(ssize_d) + strlen(ssize_i) + strlen(vect_len) + strlen(mf_n) + 5, sizeof(char));
	//+ 5 , 4 underscores, 1 for null
	strcat(file_name, ssize);
	strcat(file_name, "_");
	strcat(file_name, ssize_d);
	strcat(file_name, "_");
	strcat(file_name, ssize_i);
	strcat(file_name, "_");
	strcat(file_name, vect_len);
	strcat(file_name, "_");
	strcat(file_name, mf_n);

	return file_name;
}

/************************************************
Description:
APPENDS all epoch, trn_error,chk_error 'result.txt'.
This can be used to look at the error_result of all of epochs based on the different combination of parameters.

Inputs:
epoch_number - the current epoch
training_error - the error of training
checking_error - the error of checking

Outputs:
An output file is created in the current working directory called 'result.txt'.

Called from:
Wherever the user would like to place it.
************************************************/
void write_result(int epoch_number,int out_n, double *training_error, double *checking_error_un,double *checking_error_n, double NDEI_un, double NMSE_un, double NDEI_n, double NMSE_n, double *NMSE, double *NDEI, double *unNMSE, double *unNDEI) 
{
	int i;
	FILE *fp;

	FILE *open_file();

	fp = open_file("result.txt", "a");
	fprintf(fp, "Here are the details of final result \n");
	fprintf(fp, "The epoch number at which the training was stopped was \t %3d \n", epoch_number);
	
/*for (i=0; i< (out_n-1); i++)
	{
		fprintf(fp, "%.11f \t   %.11f\n", training_error[i], checking_error[i]);
	}*/

	//fprintf(fp, "%.11f \t   %.11f \t   %.11f \t   %.11f\n",training_error[out_n-1], checking_error[out_n-1],training_error[out_n], checking_error[out_n]);

	fprintf(fp, "At this epoch number, the training RMSE was %.11f \n", trn_rmse_error[epoch_number -1]);
fprintf(fp, "For this training RMSE & its parameters, the corresponding checking error non normalized was \t %.11f \n", checking_error_un[out_n]);
	fprintf(fp, "For this training RMSE & its parameters, the corresponding checking error normalized was \t %.11f \n", checking_error_n[out_n]);
fprintf(fp, "The corresponding NMSE non normlaized was \t %.11f \n", NMSE_un);
fprintf(fp, "The corresponding NMSE normlaized was \t %.11f \n", NMSE_n);
fprintf(fp, "The corresponding NDEI non normlized was \t %.11f \n", NDEI_un);
fprintf(fp, "The corresponding NDEI normalized was \t %.11f \n", NDEI_n);
	for(i=0; i<out_n; i++)
	{
	fprintf(fp, "for %d th variate, NMSE and NDEI normalized are %f \t %f respectively \n", i+1, NMSE[i], NDEI[i]);
	fprintf(fp, "for %d th variate, NMSE and NDEI non normalized are %f \t %f respectively \n", i+1, unNMSE[i], unNDEI[i]);
	}
	fclose(fp);
}

/**********************************************/
void calculate_root(double **data_matrix,double ** diff,int j,NODE_T **node_p) {

	double difer=0;
	int i;

	for(i = 0; i < Out_n; i++)
	{
			difer = (data_matrix[j][(i+1)*In_vect_n+i] - node_p[Node_n - (Out_n-i)]->value->real);
			diff[i][j] = difer*difer;
	}

}
/***********************************************/
void calculate_trn_err(double ** diff, double *trn_error, double *trn_datatpair_error, int training_data_n) {

	double dif=0;
	int i,j;

	for (i=0; i<Out_n; i++)
	{
		for(j = 0; j < training_data_n; j++)
		{
			dif=diff[i][j]+dif;
		}
			dif=dif/training_data_n;
			trn_error[i]=sqrt(dif);
			dif=0;
	}
	dif=0;

	for (i=0; i<training_data_n; i++)
	{
		for(j = 0; j < Out_n; j++)
		{
			dif=diff[j][i]+dif;
		}
			dif=dif/Out_n;
			trn_datapair_error[i]=sqrt(dif);
			dif=0;
	}
	dif=0;
	for (i=0; i<Out_n; i++)
	{
		for(j = 0; j < training_data_n; j++)
		{
			dif=diff[i][j]+dif;
		}
	}
	dif=dif/(training_data_n*Out_n);
	trn_error[Out_n]=sqrt(dif);

}
