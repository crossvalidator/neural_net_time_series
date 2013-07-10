#include "ancfis.h"

	int In_n, Out_n, In_vect_n, Mf_n, Rule_n, Node_n;
	int training_data_n, checking_data_n, epoch_n, parameter_n; // inputs
	int **config; //holds a binary matrix which defines node connections
	double increase_rate, decrease_rate, step_size; // inputs
	//double *step_size_pointer;
	double *trn_rmse_error, *chk_rmse_error, *trnNMSE;
	double **kalman_data;
	double *target;
	COMPLEX_T **layer_1_to_4_output;
	COMPLEX_T complexsum;
	double **training_data_matrix, **checking_data_matrix,**checking_data_matrix_un, **kalman_parameter;
	double *parameter_array, *step_size_array, **ancfis_output, **chk_output;
	double *trn_error, *chk_error_n,*chk_error_un, *trn_datapair_error, **trn_datapair_error_sorted;
	double sorting;
	int sortingindex;
	double min_trn_RMSE=100,min_chk_RMSE_n=100, min_chk_RMSE_un=100;
	int min_trn_RMSE_epoch = -1; 
	int min_chk_RMSE_epoch = -1; 
	NODE_T **node_p;
	int i,j,k,m,c;
	double tmp,tot_sqr_err,*tot_err;	//get the desired output
	FILE *fpp, *fppp, *fpppp;
	double my_array[11],*de_out,**diff;
	int bflg, aflg, errflg, ep_n;
	char *ifile;
    extern char *optarg;
    extern int optind, optopt;
	double checking_data_average_n = 0.0;
	double checking_data_average_un = 0.0;
	double checking_variance_n = 0.0;
	double checking_variance_un= 0.0;
	double NMSE_un = 0.0;
	double NMSE_n = 0.0;
	double NMSE_n2 = 0.0;
	double NDEI_un = 0.0;
	double NDEI_n = 0.0;
	double *NMSE, *NDEI, *unNMSE, *unNDEI;
	double checking_data_average_n_temp=0, temp=0, checking_data_average_un_temp =0;
	double threshold = 0.0;
	int epoch_of_threshold = -1;
	double trnvariance=0.0;
	double chkvariance = 0.0;
	double trnavg=0.0;
	double chkavg = 0.0;
	double min_trnNMSE =0.0;
#include "input.c"
#include "lib.c"
#include "complex.c"
#include "gen_config.c"
#include "datastru.c"
#include "initpara.c"
#include "debug.c"
#include "forward.c"
#include "kalman.c"
#include "de_dp.c"
#include "backward.c"
#include "chk_err.c"
#include "output.c"
#include "chaossimulatedannealling.c"
#include "new_para.c"
#include "stepsize.c"


//////////////////////////////////////////////////////////////////////////////////////////
// Define variables



/////////////////////////////////////////////////////////////////////////////////////////////
//int main(int argc, char *argv[ ])
//{
 int main(void)
 {
   	sortingindex = 0;
	if((fpp = fopen(INPUT_FILE, "r")) == NULL)
	{
		printf("Cannot open 'parameter_file'.\n");
	}

	for(j = 0; j < 11; j++)
	{
		if(fscanf(fpp, "%lf", &tmp) != EOF) 
		{
				my_array[j] = tmp;
				
		} else 		{
			printf("Not enough data in 'input_parameter'!");
		}
	}

	fclose(fpp);
 

   	In_n= (int) my_array[0];
	In_vect_n= (int) my_array[1];
	Out_n= (int) my_array[2];
	Mf_n= (int) my_array[3];
	training_data_n= (int) my_array[4];
	checking_data_n= (int) my_array[5];
	epoch_n= (int) my_array[6];
	step_size=my_array[7];
	increase_rate=my_array[8];
	decrease_rate=my_array[9];
	threshold = my_array[10];
		
	Rule_n = (int)pow((double)Mf_n, (double)In_n); //number of rules 
	Node_n = In_n + In_n*Mf_n + 3*Rule_n + In_n*Rule_n + Out_n;

	/* allocate matrices and memories */
	int trnnumcheck[training_data_n + 1];
	int trnnumchecku[training_data_n + 1];
	for(i=0; i<training_data_n +1; i++)
	{
	trnnumcheck[i]=0;
	trnnumchecku[i]=0;
	}
	
	diff =(double **)create_matrix(Out_n, training_data_n, sizeof(double)); 
	double chkvar[checking_data_n];
	double cdavg[checking_data_n];
	double chkvar_un[checking_data_n];
	double cdavg_un[checking_data_n];
	target = calloc(Out_n, sizeof(double));
	de_out = calloc(Out_n, sizeof(double));
	node_p = (NODE_T **)create_array(Node_n, sizeof(NODE_T *)); 
	config = (int **)create_matrix(Node_n, Node_n, sizeof(int)); 
	training_data_matrix = (double **)create_matrix(training_data_n, In_n*In_vect_n + Out_n, sizeof(double));
	if(checking_data_n > 0)
	{
		checking_data_matrix = (double **)create_matrix(checking_data_n, In_n*In_vect_n +Out_n, sizeof(double));
		checking_data_matrix_un = (double **)create_matrix(checking_data_n, Out_n, sizeof(double));
		chk_output =  (double **)create_matrix(checking_data_n, Out_n, sizeof(double));
	}
	layer_1_to_4_output = (COMPLEX_T **)create_matrix(training_data_n, In_n*Mf_n + 3*Rule_n, sizeof(COMPLEX_T));
	trn_rmse_error = calloc(epoch_n, sizeof(double));
	trnNMSE = calloc(epoch_n, sizeof(double));
	chk_rmse_error = calloc(epoch_n, sizeof(double));
	kalman_parameter = (double **)create_matrix(Out_n ,(In_n*In_vect_n + 1)*Rule_n, sizeof(double)); 
	kalman_data = (double **)create_matrix(Out_n ,(In_n*In_vect_n + 1)*Rule_n, sizeof(double));
	step_size_array = calloc(epoch_n, sizeof(double));
	ancfis_output = (double **)create_matrix(training_data_n , Out_n, sizeof(double)); 
	trn_error =calloc(Out_n +1, sizeof(double));
	chk_error_n = calloc(Out_n +1, sizeof(double));// changing size for adding new error measures
	chk_error_un = calloc(Out_n +1, sizeof(double));// changing size for adding new error measures
	trn_datapair_error = calloc(training_data_n, sizeof(double));
	trn_datapair_error_sorted = (double **)create_matrix(2,training_data_n, sizeof(double));
	NMSE = calloc(Out_n, sizeof(double));
	NDEI = calloc(Out_n, sizeof(double));
	unNMSE = calloc(Out_n, sizeof(double));
	unNDEI = calloc(Out_n, sizeof(double));
	//Build Matrix of 0 nd 1 to show the connected nodes
	gen_config(In_n, Mf_n,Out_n, config);//gen_config.c
	//With the above matrix, build ANCFIS connected nodes
	build_ancfis(config); //datastru.c
	//Find total number of nodes in layer 1 and 5
	parameter_n = set_parameter_mode(); //datastru.c
	parameter_array = calloc(parameter_n, sizeof(double));
	initpara(TRAIN_DATA_FILE, training_data_n, In_n, In_vect_n+1, Mf_n); // initpara.c
// after this step, the parameters (they are present in layer 1 and layer 5 only) are assigned a random initial value 
// using some basic algebra and these value are then stored in "para.ini"
	get_parameter(node_p,Node_n,INIT_PARA_FILE); //input.c
// after this step, the initial random values of the parametrs are read from "oara.ini" and assigned to the appropriate nodes in the node structure by accessing their para list.
	//Get training and testing data
	get_data(TRAIN_DATA_FILE, training_data_n, training_data_matrix); //input.c
// after this step, the training input data is read from the "data.trn" fle and stroed in the training data matrix.
	get_data(CHECK_DATA_FILE, checking_data_n, checking_data_matrix); //input.c
// after the above step, the checking data is read from the "data.chk" file and then stored in the checking data matrix.

	for(i=0; i< Out_n; i++)
	{
	for(j=0; j<training_data_n; j++)
	{
	trnavg = trnavg + training_data_matrix[j][(i+1)*In_vect_n +i];
	}
	}
	trnavg = trnavg /(Out_n * training_data_n);
	
	for(i=0; i< Out_n; i++)
	{
	for(j=0; j<training_data_n; j++)
	{
	temp = training_data_matrix[j][(i+1)*In_vect_n +i]- trnavg;
	temp = temp*temp;
	trnvariance = trnvariance + temp;
	}
	}
	trnvariance = trnvariance /((Out_n * training_data_n)-1);

	temp = 0.0;
	for(i=0; i< Out_n; i++)
	{
	for(j=0; j< checking_data_n; j++)
	{
	chkavg = chkavg + checking_data_matrix[j][(i+1)*In_vect_n +i];
	}
	}
	chkavg = chkavg /(Out_n * checking_data_n);
	
	for(i=0; i< Out_n; i++)
	{
	for(j=0; j<checking_data_n; j++)
	{
	temp = checking_data_matrix[j][(i+1)*In_vect_n +i]- chkavg;
	temp = temp*temp;
	chkvariance = chkvariance + temp;
	}
	}
	chkvariance = chkvariance /((Out_n * checking_data_n)-1);
	printf("epochs \t trn error \t tst error\n");
	printf("------ \t --------- \t ---------\n");
	//printf("not entering the epoch loop and the i loop yoyoyo\n");
/**************
	for(ep_n = 0; ep_n < epoch_n; ep_n++)
	{ 
		//step_size_pointer= &step_size;		
		//printf("epoch numbernumber %d \n", ep_n);	
		//step_size_array[ep_n] = step_size_pointer;
		step_size_array[ep_n] = step_size;
	// after the above step, the updated stepsize at the end of the last loop is stored in the step_size_array.
	// this will keep happening every time we start en epoch and hence at the end of the loop, step_size_array will 
	// have a list of all the updated step sizes. Since this is a offline version, step sizes are updated only
	// at the end of an epoch. 
		for(m = 0; m < Out_n; m++)
		{ 	
			//printf("m loop number %d \n", m);	
			for(j = 0; j < training_data_n; j++)
			{ 
				//printf("j loop number %d \n", j);				
				//copy the input vector(s) to input node(s)
				put_input_data(node_p,j, training_data_matrix); //input.c
	// after this(above) step, the input data is transferred frm the training data matrix to the "node" structure.
				//printf("testing \n");	
				//printf("reeeetesting \n");	
				target[m] = training_data_matrix[j][(m+1)*In_vect_n+m]; // *** 
	// this step assigns the value of the "m"th output of "j" th trainig data pair to target.
				//printf("testing \n");	
				//forward pass, get node outputs from layer 1 to layer 4
				calculate_output(In_n, In_n + In_n*Mf_n + 3*Rule_n - 1, j); //forward.c
	// after this step, output of nodes in layer 1 to 4 is calculated. Please note that when this happens for the first
	// time, i.e. when ep_n=0, our network parametrs are already initialized. thus, it is possible to get the
	// output of each node using the function definitios proposed in forward.c. After first epoch, our parametrs get 
	// updated and this output is then calculated using teh new parameters. The essential point to note here is that
	// we can always calculate the output of each node since we have already initialized our parameters.
				//printf("testing \n");	
				//put outputs of layer 1 to 4 into layer_1_to_4_output
		
				for(k = 0; k < Mf_n*In_n + 3*Rule_n; k++)
				{
				//printf("testing \n");	
				layer_1_to_4_output[j][k] = *node_p[k + In_n]->value;
				}
	// the above loop simply puts the values of nodes from layer 1 to layer 4 in the layer_1_to_4_output matrix.

				//identify layer 5 params using LSE (Kalman filter)
				//printf("testing \n");	
				get_kalman_data(kalman_data, target); //kalman.c
	// this function call finds out the values of O4iXnl .. these are basically the coefficients
	// of the kalman parametrs for a given training data pair
	//puts them in kalman_data matrix.
	// this kalman_data matrix has In_n number of rows and number of columns equal to number of parametrs that are
	// responsible for determining each output... as stated above, the outputs are actually the coefficients of the
	// parameters.

				//printf("testing \n");	
				//calculate Kalman parameters
				
				kalman(ep_n, j+(m*training_data_n), m, kalman_data, kalman_parameter,target); //kalman.c
	// this function call evaluates kalman parametrs for a given output, for a given epoch.. that is it takes the epoch 
	// number from us, takes the info about how many times has kalman been invoked before, also takes in the
	// output number(row number) for whihc the parametrs are to be found out... it also takes kalman_data and reads 
	// from it to estimate the kalman parameters... it also takes target .. and stores the output in the mth row of 
	// kalman_parameter.
				//printf("testing \n");	
			}
	// let me tell u what the abopve loop is doing.. after observing closely, it is easy to see that in the above loop, 
	// for a given output node, one by one, all the training data are taken inside the ANCFIS structure, outputs
	// are calculated from one to 4, then a recursive kalman filetr is used to identify the kalman
	// parametrs corresponding to the output node.. these kalman parameters are updated after every tarining data pair 
	// and finally at the end of all the training data, we have an estimate for the kalman parametrs corresponding to 		// the output node.
		}
	// thus, at the of the above loop, the kalman parametrs for all the output nodes are evaluated...

	// now, we are ready to actually calculate the outputs.. plase remember that, all this while, the actual 
	// values of the parametrs of nodes in layer 1 and layer 5 are the ones that were randomly initialized.

		for(j = 0; j < training_data_n; j++)
		{ 
			//printf("testing 1\n");	
			put_input_data(node_p,j, training_data_matrix); //input.c
			//printf("testing 2 \n");	
			for(k = 0; k < Mf_n*In_n + 3*Rule_n; k++)
			{
				*node_p[k + In_n]->value = layer_1_to_4_output[j][k];
			}
	// u must be able to see that in the above two loops, each time, whatever output we got for a given training 
	// datta pair, it was safely stored in layer_1_to_4 array...and each time, the value on the actual nodes in the
	// structure got changed.. due to new incoming training data pair..this was periodic with period trainingdata_n..
	// that is for each output node, we got the same results for a given training dat aapir.. that is the node values
	// were independent of m. Now, for a given traing data pair, we are getting those value back in the actual node 
	// node structure from that laye blh blah matrix..

			//printf("testing 3\n");	
			put_kalman_parameter(Out_n,kalman_parameter); //kalman.c
	// using this function call, we are placing the setimated value of the layer 5 parametrs in the node structure
	// by accessing each node and its parameter list.
			// Do forward pass for L5 and L6
			calculate_output(In_n + In_n*Mf_n + 3*Rule_n, Node_n, j); //forward.c
	// for a given value of the training data pair, this function calculates the output of layer 5 and layer 6 nodes 
	// and places them in the node structure.

			calculate_root(training_data_matrix,diff,j,node_p); //debug.c
	// this function call calculates the square of the erro between the predicted value of an output node and the 
	// corresponding actual value in the training data matrix..(actual output) and stores it in diff.
	// this call performs the above action for a given training data pair and for all output nodes.

			// Do backward pass and calculate all error rate for each node
			calculate_de_do(j,Out_n,target,ancfis_output); //backward.c
	// calculates de_do for each node fora given training data pair
			update_de_dp(j); //de_dp.c	
	// updates de_do for each node....
		}
	// thus at the end of this loop, estimated outputs for all the training data are calculated..also back propogatin 
	// is done and de_dp for all the nodes is updated.
		
		//printf("testing 1\n");	
		calculate_trn_err(diff,trn_error,trn_datapair_error,training_data_n); //debug.c
		//printf("testing 2 \n");	
		//training_error_measure(target,ancfis_output,diff, training_data_n, trn_error,out_n); //trn_err.c
		trn_rmse_error[ep_n] = trn_error[Out_n];
		printf("%3d \t %.11f \n", ep_n+1, trn_error[Out_n]);
		//Find RMSE of testing error
	/*************************************	if(checking_data_n != 0) 
		{
			printf("testing 3 \n");	
			epoch_checking_error(checking_data_matrix, checking_data_n, chk_error, training_data_n, chk_output, ep_n); //chk_err.c  writes to tracking.txt
			printf("testing 4 \n");	
			chk_rmse_error[ep_n] = chk_error[Out_n];
			for (i=0; i<Out_n; i++)
			//printf("%3d \t %.11f \t %.11f\n", ep_n+1, trn_error[i], chk_error[i]);
			printf("%3d \t %.11f \n", ep_n+1, trn_error[i]);
			//printf("%.11f\t %.11f\n", trn_error[Out_n],chk_error[Out_n]);
			printf("%.11f\t %.11f\n", trn_error[Out_n]);
			write_result(ep_n+1,Out_n,trn_error,chk_error);  //debug.c writes to result.txt
		} 
		else 
		{
			for (i=0; i<Out_n; i++)	
			printf("%4d \t %.11f\n", ep_n+1, trn_error[i]);
		}
***************************/

/**
		//Find minimum training error and its epoch-number
		if(trn_rmse_error[ep_n] < min_trn_RMSE) {
			min_trn_RMSE_epoch = ep_n +1;
			min_trn_RMSE = trn_rmse_error[ep_n];
			record_parameter(parameter_array);
		}

		if(ep_n < epoch_n-1)
		{ 
			//update parameters in 1st layer (Using VNCSA)
			update_parameter(1, step_size); //new_para.c
			//update stepsize
			update_step_size(trn_rmse_error, ep_n, &step_size, decrease_rate, increase_rate); //stepsize.c
		}
	}
***/
////////////////////////////////////////////////////////////

fppp = (FILE *)open_file("status.txt", "w");
fpppp = (FILE *)open_file("trn.txt", "w");


	ep_n=0;

	do
	{
		//step_size_pointer= &step_size;		
		printf("epoch numbernumber %d \n", ep_n+1);	
		//step_size_array[ep_n] = step_size_pointer;
		step_size_array[ep_n] = step_size;
	// after the above step, the updated stepsize at the end of the last loop is stored in the step_size_array.
	// this will keep happening every time we start en epoch and hence at the end of the loop, step_size_array will 
	// have a list of all the updated step sizes. Since this is a offline version, step sizes are updated only
	// at the end of an epoch. 
		for(m = 0; m < Out_n; m++)
		{ 	
			//printf("m loop number %d \n", m);	
			for(j = 0; j < training_data_n; j++)
			{ 
				//printf("j loop number %d \n", j);				
				//copy the input vector(s) to input node(s)
				put_input_data(node_p,j, training_data_matrix); //input.c
	// after this(above) step, the input data is transferred frm the training data matrix to the "node" structure.
				//printf("testing \n");	
				//printf("reeeetesting \n");	
				target[m] = training_data_matrix[j][(m+1)*In_vect_n+m]; // *** 
	// this step assigns the value of the "m"th output of "j" th trainig data pair to target.
				//printf("testing \n");	
				//forward pass, get node outputs from layer 1 to layer 4
				calculate_output(In_n, In_n + In_n*Mf_n + 3*Rule_n, j); //forward.c
	// after this step, output of nodes in layer 1 to 4 is calculated. Please note that when this happens for the first
	// time, i.e. when ep_n=0, our network parametrs are already initialized. thus, it is possible to get the
	// output of each node using the function definitios proposed in forward.c. After first epoch, our parametrs get 
	// updated and this output is then calculated using teh new parameters. The essential point to note here is that
	// we can always calculate the output of each node since we have already initialized our parameters.
				//printf("testing \n");	
				//put outputs of layer 1 to 4 into layer_1_to_4_output
		
				for(k = 0; k < Mf_n*In_n + 3*Rule_n; k++)
				{
				//printf("testing \n");	
				layer_1_to_4_output[j][k] = *node_p[k + In_n]->value;
				//fprintf(fppp, "%lf \t %lf \t \n", (layer_1_to_4_output[j][k]).real, (layer_1_to_4_output[j][k]).imag);
				}
	// the above loop simply puts the values of nodes from layer 1 to layer 4 in the layer_1_to_4_output matrix.

				//identify layer 5 params using LSE (Kalman filter)
				//printf("testing \n");	
				get_kalman_data(kalman_data, target); //kalman.c
	// this function call finds out the values of O4iXnl .. these are basically the coefficients
	// of the kalman parametrs for a given training data pair
	//puts them in kalman_data matrix.
	// this kalman_data matrix has In_n number of rows and number of columns equal to number of parametrs that are
	// responsible for determining each output... as stated above, the outputs are actually the coefficients of the
	// parameters.

				//printf("testing \n");	
				//calculate Kalman parameters
				
				kalman(ep_n, j+(m*training_data_n), m, kalman_data, kalman_parameter,target); //kalman.c
	// this function call evaluates kalman parametrs for a given output, for a given epoch.. that is it takes the epoch 
	// number from us, takes the info about how many times has kalman been invoked before, also takes in the
	// output number(row number) for whihc the parametrs are to be found out... it also takes kalman_data and reads 
	// from it to estimate the kalman parameters... it also takes target .. and stores the output in the mth row of 
	// kalman_parameter.
				//printf("testing \n");	
			}
	// let me tell u what the abopve loop is doing.. after observing closely, it is easy to see that in the above loop, 
	// for a given output node, one by one, all the training data are taken inside the ANCFIS structure, outputs
	// are calculated from one to 4, then a recursive kalman filetr is used to identify the kalman
	// parametrs corresponding to the output node.. these kalman parameters are updated after every tarining data pair 
	// and finally at the end of all the training data, we have an estimate for the kalman parametrs corresponding to 		// the output node.
		}
	// thus, at the of the above loop, the kalman parametrs for all the output nodes are evaluated...

	// now, we are ready to actually calculate the outputs.. plase remember that, all this while, the actual 
	// values of the parametrs of nodes in layer 1 and layer 5 are the ones that were randomly initialized.

		for(j = 0; j < training_data_n; j++)
		{ 
			//printf("testing 1\n");	
			put_input_data(node_p,j, training_data_matrix); //input.c
			//printf("testing 2 \n");	
			for(k = 0; k < Mf_n*In_n + 3*Rule_n; k++)
			{
				*node_p[k + In_n]->value = layer_1_to_4_output[j][k];
				/*if(ep_n==1)
				{
				fprintf(fppp, "%d.\t %lf \t + \t i%lf \n", k, (layer_1_to_4_output[j][k]).real,(layer_1_to_4_output[j][k]).imag);
				}*/
			}
	// u must be able to see that in the above two loops, each time, whatever output we got for a given training 
	// datta pair, it was safely stored in layer_1_to_4 array...and each time, the value on the actual nodes in the
	// structure got changed.. due to new incoming training data pair..this was periodic with period trainingdata_n..
	// that is for each output node, we got the same results for a given training dat aapir.. that is the node values
	// were independent of m. Now, for a given traing data pair, we are getting those value back in the actual node 
	// node structure from that laye blh blah matrix..

			//printf("testing 3\n");	
			put_kalman_parameter(Out_n,kalman_parameter); //kalman.c
			//printf("hihahahha \n");
	// using this function call, we are placing the setimated value of the layer 5 parametrs in the node structure
	// by accessing each node and its parameter list.
			// Do forward pass for L5 and L6
			calculate_output(In_n + In_n*Mf_n + 3*Rule_n, Node_n, j); //forward.c
	// for a given value of the training data pair, this function calculates the output of layer 5 and layer 6 nodes 
	// and places them in the node structure.
			//printf("hihahahha  no 2 \n");
	calculate_root(training_data_matrix,diff,j,node_p); //debug.c
	// this function call calculates the square of the erro between the predicted value of an output node and the 
	// corresponding actual value in the training data matrix..(actual output) and stores it in diff.
	// this call performs the above action for a given training data pair and for all output nodes.

			// Do backward pass and calculate all error rate for each node
			calculate_de_do(j,Out_n,target,ancfis_output); //backward.c
			//printf("hihahahha no 3 \n");
	// calculates de_do for each node fora given training data pair
			update_de_dp(j); //de_dp.c	
	// updates de_do for each node....
		}
	// thus at the end of this loop, estimated outputs for all the training data are calculated..also back propogatin 
	// is done and de_dp for all the nodes is updated.
		
		//printf("testing 1\n");	
		calculate_trn_err(diff,trn_error, trn_datapair_error, training_data_n); //debug.c
		//printf("testing 2 \n");	
		//training_error_measure(target,ancfis_output,diff, training_data_n, trn_error,out_n); //trn_err.c
		trn_rmse_error[ep_n] = trn_error[Out_n];
		trnNMSE[ep_n] = trn_rmse_error[ep_n]*trn_rmse_error[ep_n]/trnvariance;
		fprintf(fppp, "epoch number is %d \t trn RMSE is %.11f \t trn NMSE is  %lf \t \n", ep_n + 1,  trn_rmse_error[ep_n], trnNMSE[ep_n]);
		//fprintf(fpppp, "\n");
		fprintf(fpppp, "epoch number is %d \t trn RMSE is %.11f \t trn NMSE is  %lf \t \n", ep_n + 1,  trn_rmse_error[ep_n], trnNMSE[ep_n]);
		printf("trn RMSE is \t %lf \n", trn_rmse_error[ep_n]);
		printf("trn NMSE is \t %lf \n", trnNMSE[ep_n]);
		for(i=0; i<training_data_n; i++)
		{
		trn_datapair_error_sorted[0][i]=trn_datapair_error[i];
		trn_datapair_error_sorted[1][i]= i+1;
		}

		for(j=1; j<training_data_n; j++)
		{		
		for(i=0; i<training_data_n-j; i++)
		{
		if(trn_datapair_error_sorted[0][i]>trn_datapair_error_sorted[0][i+1])
		{	
		sorting=trn_datapair_error_sorted[0][i+1];
		trn_datapair_error_sorted[0][i+1]=trn_datapair_error_sorted[0][i];
		trn_datapair_error_sorted[0][i]=sorting;
		sortingindex = sorting=trn_datapair_error_sorted[1][i+1];
		trn_datapair_error_sorted[1][i+1]=trn_datapair_error_sorted[1][i];
		trn_datapair_error_sorted[1][i]=sortingindex;
		}
		}
		}

		for(j=0; j<training_data_n; j++)
		{
		fprintf(fppp, "\n");		
		fprintf(fppp, "training data pair sorted number \t %d \n", j+1);
		fprintf(fppp, "training data pair original number \t %d \n", (int)(trn_datapair_error_sorted[1][j]));
		fprintf(fppp, "training data pair sorted error in RMSE is \t %lf \n",trn_datapair_error_sorted[0][j]);
		fprintf(fpppp, "%d \t", (int)(trn_datapair_error_sorted[1][j]));
		complexsum = complex(0.0, 0.0);
		fprintf(fppp,"Normalized layer 3 outputs are as follows \n");
		for(k= In_n*Mf_n + Rule_n; k< In_n*Mf_n + 2*Rule_n; k++)
		{
		fprintf(fppp, "%d.\t %lf + i%lf \t %lf < %lf \n", k, (layer_1_to_4_output[j][k]).real,(layer_1_to_4_output[j][k]).imag, c_abs(layer_1_to_4_output[j][k]), c_phase(layer_1_to_4_output[j][k])*180/PI);
		complexsum = c_add(complexsum, layer_1_to_4_output[j][k]);
		}
		
		
		fprintf(fppp, "Sum of the outputs of layer 3 is \t %lf+i%lf \t %lf<%lf \n", complexsum.real, complexsum.imag, c_abs(complexsum), c_phase(complexsum)*180/PI);
		complexsum = complex(0.0, 0.0);
		fprintf(fppp,"dot producted layer 4 outputs are as follows \n");
		for(k=In_n*Mf_n + 2*Rule_n; k< In_n*Mf_n + 3*Rule_n; k++)
		{
		
		fprintf(fppp, "%d.\t %lf + i%lf \t %lf < %lf \n", k, (layer_1_to_4_output[j][k]).real,(layer_1_to_4_output[j][k]).imag, c_abs(layer_1_to_4_output[j][k]), c_phase(layer_1_to_4_output[j][k])*180/PI);
		complexsum = c_add(complexsum, layer_1_to_4_output[j][k]);
		}

		fprintf(fppp, "sum of the outputs of layer 4 is \t %lf +i%lf \t %lf<%lf \n", complexsum.real, complexsum.imag, c_abs(complexsum), c_phase(complexsum)*180/PI);
		if(j> training_data_n -6 )
		{
		trnnumcheck[(int)(trn_datapair_error_sorted[1][j])]= trnnumcheck[(int)(trn_datapair_error_sorted[1][j])] +1;
		}
		if(j<5 )
		{
		trnnumchecku[(int)(trn_datapair_error_sorted[1][j])]= trnnumchecku[(int)(trn_datapair_error_sorted[1][j])] +1;
		}

		}
		fprintf(fpppp, "\n");
		
		//Find RMSE of testing error
/********************************************************************************
if(checking_data_n != 0) 
		{
			printf("testing 3 \n");	
			epoch_checking_error(checking_data_matrix, checking_data_n, chk_error, training_data_n, chk_output, ep_n); //chk_err.c  writes to tracking.txt
			printf("testing 4 \n");	
			chk_rmse_error[ep_n] = chk_error[Out_n];
			for (i=0; i<Out_n; i++)
			printf("%3d \t %.11f \t %.11f\n", ep_n+1, trn_error[i], chk_error[i]);
			printf("%.11f\t %.11f\n", trn_error[Out_n],chk_error[Out_n]);
			write_result(ep_n+1,Out_n,trn_error,chk_error);  //debug.c writes to result.txt
		} 
		else 
		{
			for (i=0; i<Out_n; i++)	
			printf("%4d \t %.11f\n", ep_n+1, trn_error[i]);
		}**************************************************************************************/

		// check whether the current training RMSE is less than the threhold and store its epch number and parametrs
		
		if(trn_rmse_error[ep_n] < min_trn_RMSE) 
		{
			min_trn_RMSE_epoch = ep_n +1;
			min_trn_RMSE = trn_rmse_error[ep_n];
			min_trnNMSE = trnNMSE[ep_n];
			record_parameter(parameter_array);
		}

		if(ep_n < epoch_n-1)
		{ 
			//update parameters in 1st layer (Using VNCSA)
			update_parameter(1, step_size); //new_para.c
			//update stepsize
			update_step_size(trn_rmse_error, ep_n, &step_size, decrease_rate, increase_rate); //stepsize.c
		}
		ep_n++;
		
	} while((trnNMSE[ep_n -1]>= threshold) && (ep_n <= epoch_n -1));

for(i=1; i<=training_data_n; i++)
{
	fprintf(fpppp, "%d \t %d \n", i, trnnumcheck[i]);
}
for(i=1; i<=training_data_n; i++)
{
	fprintf(fpppp, "%d \t %d \n", i, trnnumchecku[i]);
}


if(trnNMSE[ep_n -1]< threshold)
{
fprintf(fppp, "\n");
fprintf(fppp, "We have gone below the threshold value \n");
fprintf(fppp, "the epoch number in which this happened is %d \n", min_trn_RMSE_epoch);
}
else
{
fprintf(fppp, "\n");
fprintf(fppp, "We exhausted the available epochs and threshold was not broken :( \n");
fprintf(fppp, "the epoch number which yielded minimum training RMSE is %d \n", min_trn_RMSE_epoch);
}


fclose(fppp);
fclose(fpppp);

double *minmaxc;
minmaxc= (double *)calloc(2*In_n, sizeof(double));
	
	if((fpp = fopen("minmax.txt", "r")) == NULL)
	{
		printf("Cannot open 'parameter_file'.\n");
	}

	for(j = 0; j < 2*In_n; j++)
	{
		if(fscanf(fpp, "%lf", &tmp) != EOF) 
		{
				minmaxc[j] = tmp;
				
		} else 		{
			printf("Not enough data in 'input_parameter'!");
		}
	}

	fclose(fpp);
//////////////////////////////////////////////////////////////


	restore_parameter(parameter_array); //output.c
	write_parameter(FINA_PARA_FILE); //output.c
	write_array(trnNMSE, epoch_n, TRAIN_ERR_FILE); //lib.c
	if (checking_data_n != 0)
	{
		//printf("testing 3 \n");	
		epoch_checking_error(checking_data_matrix, checking_data_n, chk_error_n, chk_error_un, training_data_n, chk_output, ep_n -1, minmaxc); //chk_err.c  writes to tracking.txt
		//printf("testing 4 \n");	
		//chk_rmse_error[ep_n] = chk_error[Out_n];
		min_chk_RMSE_n = chk_error_n[Out_n];
		printf(" initial checking RMSE is %lf \n ", min_chk_RMSE_n);
		min_chk_RMSE_un = chk_error_un[Out_n];
		//for (i=0; i<Out_n; i++)
		//printf("%3d \t %.11f \t %.11f\n", ep_n+1, trn_error[i], chk_error[i]);
		//printf("%3d \t %.11f \n", ep_n+1, trn_error[i]);
			//printf("%.11f\t %.11f\n", trn_error[Out_n],chk_error[Out_n]);
		//printf("%.11f\t \n", trn_error[Out_n]);
		//write_result(min_trn_RMSE_epoch ,Out_n,trn_rmse_error,chk_error);  //debug.c writes to result.txt about the epoch number at which the stopping was done and the corresponding training RMSE and checking RMSE
	} 
	//write_array(chk_rmse_error, epoch_n, CHECK_ERR_FILE); //lib.c
	//}
	
	write_array(step_size_array, epoch_n, STEP_SIZE_FILE); //lib.c

/**************************************************************************
	min_chk_RMSE = chk_rmse_error[epoch_n -1];
	min_chk_RMSE_epoch = epoch_n -1;	
	for(j=0; j< epoch_n; j++)
	{
	if(chk_rmse_error[j]< min_chk_RMSE)
	{
	min_chk_RMSE = chk_rmse_error[j];
	min_chk_RMSE_epoch = j;
	}
	}
*************************************************************************/
/**************************************************************
	double minmaxc[2*In_n];
	
	if((fpp = fopen("minmax.txt", "r")) == NULL)
	{
		printf("Cannot open 'parameter_file'.\n");
	}

	for(j = 0; j < 2*In_n; j++)
	{
		if(fscanf(fpp, "%lf", &tmp) != EOF) 
		{
				minmaxc[j] = tmp;
				
		} else 		{
			printf("Not enough data in 'input_parameter'!");
		}
	}

	fclose(fpp);
***************************************************************************/
	for(k=0; k< Out_n; k++)
	{
	for(j=0; j< checking_data_n; j++)
		{
		 checking_data_matrix_un[j][k]= (checking_data_matrix[j][(k+1)*In_vect_n +k])* (minmaxc[(2*k) +1] - minmaxc[2*k]) + minmaxc[2*k];
		}
	}





// the following code calculates the cdavg_un and checking datat average bothe un normalized
for(k=0; k< Out_n; k++)
	{
	for(j=0; j< checking_data_n; j++)
		{
		checking_data_average_un = checking_data_average_un + checking_data_matrix_un[j][k];
		}
	cdavg_un[k]=checking_data_average_un/checking_data_n;
	checking_data_average_un_temp=checking_data_average_un_temp+checking_data_average_un;
	checking_data_average_un=0;
		}
	
	checking_data_average_un = checking_data_average_un_temp/(Out_n*checking_data_n);
	printf("%lf is the checking datat average non normalized\n", checking_data_average_un);





// the following code calcuates the chkvar_un un normalized 
for(k=0; k< Out_n; k++)
	{	
	for(j=0; j< checking_data_n; j++)
		{				
		temp= temp + (checking_data_matrix_un[j][k] - cdavg_un[k])*(checking_data_matrix_un[j][k] - cdavg_un[k]);
		}
	chkvar_un[k]=temp/(checking_data_n-1);
	//checking_variance_un = checking_variance_un + temp;
	temp=0;
	}




temp =0.0;
// the following code cacluates the un normalized checking varinace
for(j=0; j< checking_data_n; j++)
	{
	for(k=0; k< Out_n; k++)
		{
		temp = checking_data_matrix_un[j][k] - checking_data_average_un;
		temp = temp*temp;
		checking_variance_un = checking_variance_un + temp;
		}
	}
	checking_variance_un = checking_variance_un/((Out_n*checking_data_n)-1);
	printf("%lf is the checking variance non normalized \n", checking_variance_un);

temp =0.0;




checking_data_average_n=0.0;
checking_data_average_n_temp=0.0;
// the following code calculates the cdavg and checking data average both normalized
for(k=0; k< Out_n; k++)
	{	
for(j=0; j< checking_data_n; j++)
		{		
		checking_data_average_n = checking_data_average_n + checking_data_matrix[j][(k+1)*In_vect_n +k];
		}
		cdavg[k]=checking_data_average_n/checking_data_n;
		checking_data_average_n_temp=checking_data_average_n_temp+checking_data_average_n;
		checking_data_average_n=0;
	}
	checking_data_average_n = checking_data_average_n_temp/(Out_n*checking_data_n);
	printf("%lf is the checking datat average  normalized\n", checking_data_average_n);


temp =0.0;
checking_variance_n =0.0;
// the following code cacluates the normalized checking varinace
for(j=0; j< checking_data_n; j++)
	{
	for(k=0; k< Out_n; k++)
		{
		temp = checking_data_matrix[j][(k+1)*In_vect_n +k] - checking_data_average_n;
		temp = temp*temp;
		checking_variance_n = checking_variance_n + temp;
		}
	}
checking_variance_n = checking_variance_n/((Out_n*checking_data_n)-1);
temp = 0.0;
printf("%lf is the checking variance normalized \n", checking_variance_n);



// the following code calcuatres the normalized chkvar[k]
temp=0.0;
for(k=0; k< Out_n; k++)
	{
	for(j=0; j< checking_data_n; j++)
		{	
		temp= temp + (checking_data_matrix[j][(k+1)*In_vect_n +k] - cdavg[k])*(checking_data_matrix[j][(k+1)*In_vect_n +k] - cdavg[k]);
	}
	chkvar[k]=temp/(checking_data_n-1);
	//checking_variance_n = checking_variance_n + temp;
	temp=0;
	}


	
	
	

	NMSE_un = min_chk_RMSE_un * min_chk_RMSE_un / checking_variance_un;
	NMSE_n = min_chk_RMSE_n * min_chk_RMSE_n / checking_variance_n;
	NMSE_n2 = min_chk_RMSE_n * min_chk_RMSE_n / chkvariance;
	NDEI_un = sqrt(NMSE_un);
	NDEI_n = sqrt(NMSE_n);




	for(k=0;k<Out_n;k++)
	{
	NMSE[k]=chk_error_n[k]*chk_error_n[k]/chkvar[k];
	NDEI[k]=sqrt(NMSE[k]);
	unNMSE[k]=chk_error_un[k]*chk_error_un[k]/chkvar_un[k];
	unNDEI[k]=sqrt(unNMSE[k]);
	}






	write_result(min_trn_RMSE_epoch ,Out_n,trn_rmse_error,chk_error_un, chk_error_n, NDEI_un, NMSE_un, NDEI_n, NMSE_n, NMSE, NDEI, unNMSE, unNDEI); //debug.c writes to result.txt about the epoch number at which the stopping was done and the corresponding training RMSE and checking RMSE
	printf("Minimum training RMSE is \t %f \t \n",min_trn_RMSE); 
	printf("Minimum training RMSE epoch is \t %d \n",min_trn_RMSE_epoch); 
	printf("Minimum training NMSE is \t %f \t \n",min_trnNMSE); 
	//printf("Minimum training RMSE epoch is \t %d \n",min_trnNMSE_epoch); 
	//printf("Minimum training RMSE is \t %f \t \n",min_trn_RMSE); 
	//printf("Minimum training RMSE epoch is \t %d \n",min_trn_RMSE_epoch); 
	printf("%f \t is the checking RMSE non normalized\n",min_chk_RMSE_un);
	printf("%f \t is the checking RMSE normalized\n",min_chk_RMSE_n);
	//printf("%f \t is the checking RMSE normalized22222222 \n",min_chk_RMSE_n2);
	printf(" checking NMSE non normlized is %f \t NDEI non normalized is %f \n",NMSE_un, NDEI_un); 
	printf("checking NMSE normalized is %f \t NDEI normalized is %f \n",NMSE_n, NDEI_n); 
	printf("checking NMSE2 normalized is %f \n",NMSE_n2); 
	printf("traning data variance is  %f \n",trnvariance); 
	return(0);
}
