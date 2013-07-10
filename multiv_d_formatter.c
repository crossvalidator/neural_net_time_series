/*
This is a stand alone program takes in a time series data file and formats it for use in ANCFIS.
The input data file should gave only 1 column, each entry is on its own row.

This should be compiled and called from Matlab.

To compile:
type "mex data_formatter.c" without the quotation marks.

An example of calling the program:
n_trn_data = data_formatter('sunspot.txt','data.trn',305, 120, 1, 1, 1, 1, 1, 1, 0);
A description of each input can be seen in the comments below.

For an example of what some of the parameters mean,
suppose we have a file:
1
2
3
4
5
6
7
8
9
10

ip_vector_length = 4 with all else at 1 will give:
1 2 3 4|5
2 3 4 5|6
3 4 5 6|7
4 5 6 7|8
5 6 7 8|9
6 7 8 9|10

(note that the | indicates a seperation between inputs and outputs,
it does not actually appear in the file)

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

#include "mex.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#define MIN_REAL      -HUGE_VAL
#define MAX_REAL      +HUGE_VAL
#define LO            0.1
#define HI            0.9
#define MIN(x,y)      ((x)<(y) ? (x) : (y))
#define MAX(x,y)      ((x)>(y) ? (x) : (y))

/***********************************************
Description:
This is a seperate program that is used in conjunction with ancfismex.
This program takes in a single columned time series data file and formats it
for proper use in ANCFIS.

Generally, this program MUST be used before trying to run ANCFIS on a data file.

Inputs:
nrhs - Should be 11.
prhs[0] - The desired input file name.
	The input file must be a SINGLE column or 1 entry per line.
prhs[1] - The desired name for the output file.
	For ANCFIS, this should be 'data.trn" or "data.chk".  These defaults can be changed in ancfis.h.
prhs[2] - The number of rows in the input file.
prhs[3] - The desired input vector length.
prhs[4] - The desired space between the elements of the input vector.
	This can be used for a lagged representation of the input vector.
prhs[5] - The space between the last input and the first output.
	This can be used to specify a delay in predicting outputs.
prhs[6] - The length of the output vector.
prhs[7] - The spacing between the elements of the output vector.
prhs[8] - The space between the vectors themselves.
prhs[9] - 1 means normalize the inputs, 0 means do not normalize the inputs
prhs[10] - 1 means normalize the outputs, 0 means do not normalize the outputs

Outputs:
nlhs - Should be 1.
plhs[0] - Returns the nubmer of lines (or data pairs) in the output file.
	This is needed since ancfismex needs to know the number of lines in the training/checking data files.

Called from:
Matlab
************************************************/
void
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
const mxArray *prhs[];
mxArray *plhs[];
{
	int i, j, total_span, num_iterations;
	char *first_file_name, *second_file_name,*target_file_name,*op_file_name;
  	int buflen, status, number_of_rows, ip_vector_length, ip_spacing, ip_op_spacing, op_vector_length, op_spacing, pair_spacing, ip_normalize, op_normalize;
	FILE *first_input_file, *second_input_file, *target_input_file,*output_file;
	double *first_trn_data, *second_trn_data, *target_trn_data, max, *matlab_output,min,tmp;

    if (nrhs != 13) {
	     mexErrMsgTxt("Usage: first_file_name, second_file_name, target_file_name, op_file_name, number_of_rows_in_file, input_vector_length, input_spacing, space_between_last_input_and_first_output, output_vector_length, output_spacing, space_between_data_pairs, ip_normalize, op_normalize\n");
	}

	buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
	first_file_name = mxCalloc(buflen, sizeof(char));

	status = mxGetString(prhs[0], first_file_name, buflen);
	if (status != 0) {
    	mexWarnMsgTxt("Not enough space. String is truncated.");
	}

    buflen = (mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;
	second_file_name = mxCalloc(buflen, sizeof(char));

	status = mxGetString(prhs[1], second_file_name, buflen);
	if (status != 0) {
    	mexWarnMsgTxt("Not enough space. String is truncated.");
	}

    buflen = (mxGetM(prhs[2]) * mxGetN(prhs[2])) + 1;
	target_file_name = mxCalloc(buflen, sizeof(char));

	status = mxGetString(prhs[2], target_file_name, buflen);
	if (status != 0) {
    	mexWarnMsgTxt("Not enough space. String is truncated.");
	}

	buflen = (mxGetM(prhs[3]) * mxGetN(prhs[3])) + 1;
	op_file_name = mxCalloc(buflen, sizeof(char));

	status = mxGetString(prhs[3], op_file_name, buflen);
	if (status != 0) {
    	mexWarnMsgTxt("Not enough space. String is truncated.");
	}

	number_of_rows = (int) mxGetScalar(prhs[4]); //min 1
	ip_vector_length = (int) mxGetScalar(prhs[5]); //min 1
	ip_spacing = (int) mxGetScalar(prhs[6]); //min 1
	ip_op_spacing = (int) mxGetScalar(prhs[7]); //min 1
	op_vector_length = (int) mxGetScalar(prhs[8]); //min 1
	op_spacing = (int) mxGetScalar(prhs[9]); //note this is irrelevant if op_vector_length is 1
	pair_spacing = (int) mxGetScalar(prhs[10]); //min 1
	ip_normalize = (int) mxGetScalar(prhs[11]); //0 or 1
	op_normalize = (int) mxGetScalar(prhs[12]); //0 or 1

	if (number_of_rows < 1 || ip_vector_length < 1 || ip_spacing < 1 || ip_op_spacing < 1 || op_vector_length < 1 || op_spacing < 1 || pair_spacing < 1) {
		mexErrMsgTxt("Spacing parameters cannot be < 1.\n");
	}

	if ((first_input_file = fopen(first_file_name, "r")) == NULL) {
		mexErrMsgTxt("Cannot open 'first_data_file'.\n");
	}
	if ((second_input_file = fopen(second_file_name, "r")) == NULL) {
		mexErrMsgTxt("Cannot open 'second_data_file'.\n");
	}
    if ((target_input_file = fopen(target_file_name, "r")) == NULL) {
		mexErrMsgTxt("Cannot open 'target_data_file'.\n");
	}
	mxFree(first_file_name);
		mxFree(second_file_name);
			mxFree(target_file_name);

	min = MAX_REAL;
	max = MIN_REAL;


	//store the time series data in an array
	first_trn_data = (double *)mxCalloc(number_of_rows, sizeof(double));
	for (i = 0 ; i < number_of_rows ; i++) {
		fscanf(first_input_file, "%lf", &tmp);
		first_trn_data[i] = tmp;
 
		min = MIN(min, first_trn_data[i]);
		max = MAX(max, first_trn_data[i]);
	}

	fclose(first_input_file);

	//mexPrintf("%f, %f", min, max);
	if (ip_normalize == 1) {
		if (min < 0) {
			for (i = 0 ; i < number_of_rows ; i++) {
				first_trn_data[i] += fabs(min);
			}

			max += fabs(min);
		}
			
		if(max != min) {
			for (i = 0 ; i < number_of_rows ; i++) {
				//first_trn_data[i] = (first_trn_data[i]-min)/(max-min)* (HI-LO) + LO;
				first_trn_data[i] /= max;
			}
		} else {
			mexErrMsgTxt("ERR in normalize");
		}
	}

	min = MAX_REAL;
	max = MIN_REAL;


	second_trn_data = (double *)mxCalloc(number_of_rows, sizeof(double));
	for (i = 0 ; i < number_of_rows ; i++) {
		fscanf(second_input_file, "%lf", &tmp);
		second_trn_data[i] = tmp;
 
    	min = MIN(min, second_trn_data[i]);
		max = MAX(max, second_trn_data[i]);
	}

	fclose(second_input_file);


	if (ip_normalize == 1) {
		if (min < 0) {
			for (i = 0 ; i < number_of_rows ; i++) {
				second_trn_data[i] += fabs(min);
			}

			max += fabs(min);
		}		
	
		if(max != min) {
			for (i = 0 ; i < number_of_rows ; i++) {
				//second_trn_data[i] = (second_trn_data[i]-min)/(max-min)* (HI-LO) + LO;
				second_trn_data[i] /= max;
			}
		} else {
			mexErrMsgTxt("ERR in normalize");
		}
	}


	min = MAX_REAL;
	max = MIN_REAL;


	target_trn_data = (double *)mxCalloc(number_of_rows, sizeof(double));
	for (i = 0 ; i < number_of_rows ; i++) {
		fscanf(target_input_file, "%lf", &tmp);
		target_trn_data[i] = tmp;
 
    		min = MIN(min, target_trn_data[i]);
	    	max = MAX(max, target_trn_data[i]);
	}

	fclose(target_input_file);


	if (ip_normalize == 1) {
		if (min < 0) {
			for (i = 0 ; i < number_of_rows ; i++) {
				target_trn_data[i] += fabs(min);
			}

			max += fabs(min);
		}		
	
		if(max != min) {
			for (i = 0 ; i < number_of_rows ; i++) {
				//target_trn_data[i] = (target_trn_data[i]-min)/(max-min)* (HI-LO) + LO;
				target_trn_data[i] /= max;
			}
		} else {
			mexErrMsgTxt("ERR in normalize");
		}
	}

	//calculate the number of iterations

	//need to get max span of 1 iteration
	//span of input is (ip_vector_length-1)*ip_spacing
	//span b/w input and output is ip_op_spacing
	//span b/w outputs is (op_vector_length-1)*op_spacing
	//total_span = those 3 spans + 1

	//num_iterations = 1 + (number_of_rows-total_span)/pair_spacing
	total_span = 1 + (ip_vector_length - 1)*ip_spacing + ip_op_spacing + (op_vector_length - 1)*op_spacing;
	num_iterations = 1 + floor((number_of_rows - total_span)/pair_spacing);

	if ((output_file = fopen(op_file_name, "w")) == NULL) {
		mexErrMsgTxt("Cannot write training data.\n");
	}

	mxFree(op_file_name);

	for (i = 0 ; i < num_iterations ; i++) {
		//the input vector
		for (j = 0 ; j < ip_vector_length ; j++) {
			fprintf(output_file, "%4.10lf ", first_trn_data[i*pair_spacing+j*ip_spacing]);
		}
		
		for (j = 0 ; j < ip_vector_length ; j++) {
			fprintf(output_file, "%4.10lf ", second_trn_data[i*pair_spacing+j*ip_spacing]);

		}
		
		//the output vector
		for(j = 0 ; j < op_vector_length ; j++) {
			if(op_normalize == 0 && ip_normalize == 1) { //undo the normalization
				tmp = target_trn_data[i*pair_spacing+(ip_vector_length - 1)*ip_spacing + ip_op_spacing + j*op_spacing]*(max-min)+min;
			} else {
				tmp = target_trn_data[i*pair_spacing+(ip_vector_length - 1)*ip_spacing + ip_op_spacing + j*op_spacing];
			}

			if(j != op_vector_length - 1) {
				fprintf(output_file, "%4.10lf ", tmp);
			} else {
				fprintf(output_file, "%4.10lf", tmp);
			}
		}

		if(i != num_iterations - 1) {
			fprintf(output_file, "\n");
		}
	}

	fclose(output_file);
	mxFree(first_trn_data);
	mxFree(second_trn_data);
	mxFree(target_trn_data);

	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	matlab_output = mxGetPr(plhs[0]);
	*matlab_output = num_iterations;
}
