//#include "ancfis.h"
/*
This file handles the code for the adaptive step size.

There are 3 functions:
1. update_step_size
	2. check_increase_ss
	3. check_decrease_ss
*/

/***********************************************
Description:
Update the step size.

Inputs:
*error_array - array with previous RMSE training errors
i - epoch number
*step_size_p - pointer to a double value step size (not an array)
decrease_rate - the stepsize decrease rate
increase_rate - the stepsize increase rate

Outputs:
Directly updates *step_size.

Called from:
mexFunction in ancficmex.c
************************************************/
void update_step_size(double *error_array, int i, double *step_size_p, double decrease_rate, double increase_rate) {
	static int last_decrease_ss, last_increase_ss;
	int check_increase_ss();
	int check_decrease_ss();

	if(i == 0)
		last_decrease_ss = last_increase_ss = 0;

	if(check_decrease_ss(error_array, last_decrease_ss, i)) {
		*step_size_p *= decrease_rate;
		//printf("Ss decrease to %f after epoch %d.\n", *step_size_p, i+1);
		last_decrease_ss = i;
		return;
	}

	if(check_increase_ss(error_array, last_increase_ss, i)) {
		*step_size_p *= increase_rate;
		//printf("Ss increase to %f after epoch %d.\n", *step_size_p, i+1);
		last_increase_ss = i;
		return;
	}
}

/***********************************************
Description:
Check whether we should increase the step size.  The rule for doing this is
4 consecutive decreases in the training RMSE measure.

Inputs:
*error_array - array with previous RMSE training errors
last_change - stores the epoch of the previous increase
current - current epoch

Outputs:
Return 1 if the step size needs to be increased, 0 otherwise.

Called from:
update_step_size in stepsize.c
************************************************/
int check_increase_ss(double *error_array, int last_change, int current) {
	if(current - last_change < 4)
		return(0);
	if((error_array[current] < error_array[current - 1]) &&
	    (error_array[current - 1] < error_array[current - 2]) &&
	    (error_array[current - 2] < error_array[current - 3]) &&
	    (error_array[current - 3] < error_array[current - 4]))
		return(1);
	return(0);
}

/***********************************************
Description:
Check whether we should decrease the step size.  The rule for doing this is
4 consecutive alternating RMSE measures.

Inputs:
*error_array - array with previous RMSE training errors
last_change - stores the epoch of the previous decrease
current - current epoch

Outputs:
Return 1 if the step size needs to be decreased, 0 otherwise.

Called from:
update_step_size in stepsize.c
************************************************/
int check_decrease_ss(double *error_array, int last_change, int current) {
	if(current - last_change < 4)
		return(0);
	if((error_array[current] < error_array[current - 1]) &&
	    (error_array[current - 1] > error_array[current - 2]) &&
	    (error_array[current - 2] < error_array[current - 3]) &&
	    (error_array[current - 3] > error_array[current - 4]))
		return(1);
	return(0);
}

