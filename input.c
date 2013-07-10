
//#include "ancfis.h"
/*
int In_n, Out_n, In_vect_n, Mf_n, Rule_n, Node_n;
	int training_data_n, checking_data_n, epoch_n, parameter_n; // inputs
	int **config; //holds a binary matrix which defines node connections
	double increase_rate, decrease_rate, step_size; // inputs
	double *trn_rmse_error, *chk_rmse_error;
	double **kalman_data;
	double *target;
	COMPLEX_T **layer_1_to_4_output;
	double **training_data_matrix, **checking_data_matrix, **kalman_parameter;
	double *parameter_array, *step_size_array, **ancfis_output, **chk_output;
	double *trn_error, *chk_error;
	double min_trn_RMSE=100,min_chk_RMSE=100;
	double min_trn_RMSE_epoch = -1; 
	double min_chk_RMSE_epoch = -1; 
	NODE_T **node_p;
	int i,j,k,m,c;
	double tmp,tot_sqr_err,*tot_err;	//get the desired output
	FILE *fpp;
	double my_array[10],*de_out,**diff;
	int bflg, aflg, errflg, ep_n;
	char *ifile;
    extern char *optarg;
    extern int optind, optopt;

*/













/***********************************************
Description:
Get the initial parameters for layer 1 and 4 from "para.ini" and place them
into the ANCFIS structure.
"para.ini" is generated from initpara.c

Inputs:
*parameter_file - the name of the file "para.ini"

Outputs:
Operates on global variable node_p

Called from:
mexFunction in ancfismex.c
************************************************/
void get_parameter(NODE_T **node_p,int Node_n,char *parameter_file)
 {
	int i, j;
	int parameter_n;
	double tmp;
	FILE *fp;
	PARAMETER_LIST_T *p;

	if((fp = fopen(parameter_file, "r")) == NULL)
	{
		printf("Cannot open 'parameter_file'.\n");
	}
//	printf(" %d:",Node_n);

	for(i = 0; i < Node_n; i++) {
		if(node_p[i]->parameter == NULL)
			continue;

 	parameter_n = node_p[i]->parameter_n;
	p = node_p[i]->parameter;


		for(j = 0; j < parameter_n; j++) {
			if(fscanf(fp, "%lf", &tmp) != EOF) {
				p->content = tmp;
				p = p->next;
			} else {
				printf("Not enough data in 'init_parameter'!");
			}
		}
	}

	fclose(fp);
}

/***********************************************
Description:
Transfer (training or checking) data from data_file to matrix 'data_matrix'.
This is needed for an "internal representation" of the data file instead of having to
open up the file every time we need to access data.

Inputs:
*data_file - the name of the training of checking data text file (from data_formatter.c)
data_n - the number of rows in the data file
**data_matrix - a pre initialized matrix that will hold the training or checking data

Outputs:
data_matrix is filled with the training or checking data

Called from:
mexFunction in ancfismex.c
************************************************/
void get_data(char *data_file, int data_n, double **data_matrix) 
{
	int i, j;
	double tmp;
	FILE *fp;

	if(data_n == 0)
		return;

	if((fp = fopen(data_file, "r")) == NULL)
	{
		printf("Cannot open 'data_file'.\n");
	}

	for(i = 0; i < data_n; i++)
	{
		for(j = 0; j < (In_n*In_vect_n) +Out_n ; j++) {
			if(fscanf(fp, "%lf", &tmp) != EOF) {
				data_matrix[i][j] = tmp;
			} else {
				printf("Not enough data!");
			}
		}
	}
	fclose(fp);
}

/***********************************************
Description:
Put input part of (j+1)-th data (training or checking) to input nodes.
Operates on global variable node_p

Inputs:
j - the current index of the training or checking data pair
**data_matrix - training or checking data matrix

Outputs:
node_p->value = the norm of the data window
node_p->trn_data = the data window

Called from:
mexFunction in ancfismex.c
epoch_checking_error in chk_err.c

************************************************/
void put_input_data(NODE_T **node_p,int j, double **data_matrix)
{
	int i, k;
	double *input_vect;
	input_vect = calloc(In_vect_n, sizeof(double));

	COMPLEX_T c=complex(0.0,0.0);

	for(i = 0 ; i < In_n ; i++) 
	{
		//input_vect = calloc(In_vect_n, sizeof(double));

		for(k = 0 ; k < (In_vect_n) ; k++)
		{
			input_vect[k] = data_matrix[j][k + (i*(In_vect_n+1))];
		}

		*node_p[i]->value = c;
		node_p[i]->input_vector = input_vect;
		//free(input_vect);
	}

	free(input_vect);

}
