
// please read only this instruction.
// this is new code after discussion with prof scott dick. here, layer 3 and layer 4 are redefined.
// this is improvemen no 2


// this is a basic code of teh original entire code... many other improvements have been made in this code and these 
// improved forward.c codes are typed below this code with proper commenting. Please have a look at them too.
//#include "ancfis.h"
//#include "complex.c"
/*
This file contains the functions which calculate the outputs for the forward pass.
In our case, there are 6 layers.  Each layer has its own function.

Layer 1 -> mf
Layer 2 -> multiply
Layer 3 -> normalize
Layer 4 -> dot product
Layer 5 -> consequent
Layer 6 -> sum

If any of these functions are changed, make sure the corresponding
functions in backward.c or de_dp.c are changed as well.

There are 8 functions:
1. calculate_output
	2. input
	3. mf
	4. multiply
	5. normalize
	6. dotproduct
	7. consequent
	8. sum
*/
/***********************************************
Description:
Calculate node outputs from node 'from' to node 'to' for the forward pass.
Using this function, we can calculate the output for a single node or
for whole layers.

This function will automatically call the 5 forward pass functions.
The individual forward pass functions should not be called on their own.

Inputs:
from - starting node index
to - ending node index
iteration_number - current training data pair number

Outputs:
Operates on the global variable 'node_p'

Called from:
mexFunction in ancfismex.c
************************************************/

void calculate_output(int from, int to, int iteration_number) {
	COMPLEX_T input(), mf(), multiply();
	COMPLEX_T normalize(), dotproduct(), consequent(), sum();
	static COMPLEX_T (*function[7])() = {input, mf, multiply, normalize, dotproduct, consequent, sum};
	int i,j;
	int function_index;

	if((from > to) || (from < In_n) || (to > Node_n))
		printf("Error in calculate_output!");

	for(i = from ; i < to ; i++) {
		function_index = node_p[i]->function_index;

		*node_p[i]->value = (*function[function_index])(i, iteration_number);

			
	}

}

/***********************************************
Description:
This is a placeholder function that should never be called.

Inputs:
N/A

Outputs:
N/A

Called from:
N/A
************************************************/

struct complex_s input(int node_index, int iteration_number) {
	printf("This shouldn't have been called!\n");
	return(complex(0.0, 0.0));
}

/***********************************************
Description:
Calculates the output of layer 1.
This is currently defined as the sum of the convolution
between the data set and the sampled sine membership function.
Note that the convolution isn't actually calculated since
we are taking the summation at the end.  This is equivalent to just
summing every single combination of the 2 data sets.

The convolution sum is then applied to the Elliot function.

Inputs:
node_index - the node index
iteration_number - the current training data pair number

Outputs:
The sum of the convolution

Called from:
calculate_output in forward.c
************************************************/

COMPLEX_T mf(int node_index, int iteration_number) {
	NODE_LIST_T *arg_p = node_p[node_index]->fan_in;
	PARAMETER_LIST_T *para_p = node_p[node_index]->parameter;
	double a, b, c, d, t, r, theta;
	double *trn_data2 = arg_p->content->input_vector;
	int i, j;
	COMPLEX_T sum_conv, tmp, tmp2;

	a = para_p->content;
	b = para_p->next->content;
	c = para_p->next->next->content;
	d = para_p->next->next->next->content;
	t = 2*PI/In_vect_n;

	sum_conv = complex(0.0, 0.0);
	
	for(i = 0 ; i < In_vect_n ; i++) {
		for(j = In_vect_n - 1 ; j >= 0 ; j--) { //going backwards not needed either
			theta = t*(j + iteration_number);
			r = d*sin(a*theta + b) + c;
			tmp = complex(r*cos(theta), r*sin(theta));
			tmp2 = c_mul_scalar(tmp, trn_data2[i]);
			sum_conv = c_add(sum_conv, tmp2);
		}
	}

    tmp = c_mul_scalar(sum_conv, 1/(c_abs(sum_conv) + 1));
	return(tmp);
}

/***********************************************
Description:
Calculates the output of layer 2.
Multiplies the output values of all the fan-in nodes
to the current node index.

Inputs:
node_index - the node index
iteration_number - the current training data pair number (not used)

Outputs:
The product of all the fan-in nodes to this one.

Called from:
calculate_output in forward.c
************************************************/

struct complex_s multiply(int node_index, int iteration_number) {
	NODE_LIST_T *p;
	NODE_LIST_T *arg_p = node_p[node_index]->fan_in;
	PARAMETER_LIST_T *para_p = node_p[node_index]->parameter;
	COMPLEX_T product = complex(1.0, 0.0);

	if(para_p != NULL)
		printf("Error in multiply!");   // we should have actually checked for layer of the node_index.. that would have been more confirmative
		
	for(p = arg_p ; p != NULL ; p = p->next) {
		product = c_mul(product, *p->content->value);
	}

	return(product);
}

/***********************************************
Description:
Calculates the output of layer 3.
Normalizes the output based on the values of layer 2.

Inputs:
node_index - the node index
iteration_number - the current training data pair number (not used)

Outputs:
The normalized output.

Called from:
calculate_output in forward.c
************************************************/

COMPLEX_T normalize(int node_index, int iteration_number) {
	NODE_LIST_T *arg_p = node_p[node_index]->fan_in;
	PARAMETER_LIST_T *para_p = node_p[node_index]->parameter;
	int i;
	COMPLEX_T tmp;
	NODE_LIST_T *p;
	COMPLEX_T actual_sum = complex(0.0, 0.0);

	if(para_p != NULL)
		printf("Error in normalize!");

	for(p = arg_p ; p != NULL ; p = p->next) { //loop for all fan in
		actual_sum = c_add(actual_sum, *p->content->value);
	}

	p = arg_p;
	for(i = 0 ; i < node_p[node_index]->local_index ; i++) {  //shouldn't this be local_index + 1 .. will have to check the definition
		p = p->next;
	} //after this loop p is the current weight wi

    tmp = c_div(*p->content->value, actual_sum);
	//printf("normalize \n");
	//printf("%lf \t %lf \t \n", tmp.real, tmp.imag);
    return (tmp);
}

/***********************************************
Description:
Calculates the dot product of the current weight wi
with the sum of all previous weights.
The result is always real valued, but we return a complex number
(with the imaginary portion set to zero).

Inputs:
node_index - the node index
iteration_number - the current training data pair number (not used)

Outputs:
The normalized output.

Called from:
calculate_output in forward.c
************************************************/
/***********************************************
COMPLEX_T dotproduct(int node_index, int iteration_number) {
	NODE_LIST_T *arg_p = node_p[node_index]->fan_in;
	PARAMETER_LIST_T *para_p = node_p[node_index]->parameter;
	int i;
	//COMPLEX_T weight_sum = complex(0.0, 0.0);
	NODE_LIST_T *p;
	COMPLEX_T dum;

	//printf("dot product \n");
	if(para_p != NULL)
		printf("Error in normalize1!");
	/*
	for(p = arg_p ; p != NULL ; p = p->next) {
		weight_sum = c_add(weight_sum, *p->content->value);
	}
	*/
	//printf("%lf \t %lf \t \n", weight_sum.real, weight_sum.imag);
/*
    p = arg_p;
	for(i = 0 ; i < node_p[node_index]->local_index ; i++) {
		p = p->next;
	} //after this loop p is the current weight wi
	//printf("%lf \t %lf \t \n", (*p->content->value).real, (*p->content->value).imag);
	dum = complex((*p->content->value).real, 0.0);
	
	//printf("%lf \t %lf \t \n", dum.real, dum.imag);
    return (dum);
}
********************************************/
/****************************************
COMPLEX_T dotproduct(int node_index, int iteration_number) {
	NODE_LIST_T *arg_p = node_p[node_index]->fan_in;
	PARAMETER_LIST_T *para_p = node_p[node_index]->parameter;
	int i;
	COMPLEX_T weight_sum = complex(0.0, 0.0);
	NODE_LIST_T *p;
	COMPLEX_T dum;

	//printf("dot product \n");
	if(para_p != NULL)
		printf("Error in normalize1!");

	for(p = arg_p ; p != NULL ; p = p->next) {
		weight_sum = c_add(weight_sum, *p->content->value);
	}
	//printf("%lf \t %lf \t \n", weight_sum.real, weight_sum.imag);
    p = arg_p;
	for(i = 0 ; i < node_p[node_index]->local_index ; i++) {
		p = p->next;
	} //after this loop p is the current weight wi
	//printf("%lf \t %lf \t \n", (*p->content->value).real, (*p->content->value).imag);
	dum = c_dot_product(*p->content->value, weight_sum);
	
	//printf("%lf \t %lf \t \n", dum.real, dum.imag);
    return (dum);
}
*********************************************/

COMPLEX_T dotproduct(int node_index, int iteration_number) {
	NODE_LIST_T *arg_p = node_p[node_index]->fan_in;
	PARAMETER_LIST_T *para_p = node_p[node_index]->parameter;
	int i;
	
	NODE_LIST_T *p;
	double mag_sum = 0.0;
        //COMPLEX_T gen_sum = complex(0.0, 0.0);
      //  double abs_weight_sum = 0.0;
	//double sqr_mag_sum = 0.0;
	//double sqrof_abs_weight_sum = 0.0;
	COMPLEX_T weight_sum = complex(0.0, 0.0);

	if(para_p != NULL)
		printf("Error in normalize1!");

	for(p = arg_p ; p != NULL ; p = p->next) {
		weight_sum = c_add(weight_sum, *p->content->value);
		//sqr_mag_sum = sqr_mag_sum + (c_abs(*p->content->value)*c_abs(*p->content->value));
		
	}
	 //abs_weight_sum = c_abs(weight_sum);
	  //sqrof_abs_weight_sum = abs_weight_sum * abs_weight_sum;
	
    p = arg_p;
	for(i = 0 ; i < node_p[node_index]->local_index ; i++) {
		p = p->next;
	} //after this loop p is the current weight wi



return (c_dot_product(*p->content->value, c_sub(weight_sum,*p->content->value)));
}


/***********************************************
Description:
Calculates the output of layer 5.
Note that the Kalman filter algorithm must be run to identify
the consequent parameters {p, r} before calling this.  {p, r} are
currently real values.

Inputs:
node_index - the node index
iteration_number - the current training data pair number (not used)

Outputs:
w(px+r), for 1 input.  w and x are complex but with the imaginary part set to 0.

Called from:
calculate_output in forward.c
************************************************/

COMPLEX_T consequent(int node_index, int iteration_number) {
	NODE_LIST_T *arg_p = node_p[node_index]->fan_in;
	PARAMETER_LIST_T *para_p = node_p[node_index]->parameter;
	int i, h;
	double a, tmp, total=0;
	COMPLEX_T x;

	for(i = 0; i < In_n + 1; i++) {
		double *trn_data = node_p[i]->input_vector;
	    x = *arg_p->content->value;
		//x = 0 for the fan in nodes of the input layer
		//x = wi on last iteration, output of layer 4

		if(i == In_n)
			break;

		for(h = 0 ; h < In_vect_n ; h++) {
			a = para_p->content;//this is p, or q
			tmp = a*(trn_data[h]); //eg p1x1, p2x2, q1x1, q2x2
			total+=tmp;
			para_p = para_p->next;
		}
		arg_p = arg_p->next; //get next fan in node from layer 4
	}
	a = para_p->content; //this is r
    tmp = total+a;
	return(c_mul_scalar(x, tmp));// 
}

/***********************************************
Description:
Calculates the output of layer 6.
Sums all the outputs of the fan in nodes to this node.

Inputs:
node_index - the node index
iteration_number - the current training data pair number (not used)

Outputs:
Final predicted value.

Called from:
calculate_output in forward.c
************************************************/

COMPLEX_T sum(int node_index, int iteration_number) {
	NODE_LIST_T *arg_p = node_p[node_index]->fan_in;
	NODE_LIST_T *t;
	COMPLEX_T total = complex(0.0, 0.0);

	if(arg_p == NULL)
		printf("Error! Given pointer is NIL!");

	for(t = arg_p; t != NULL; t = t->next) {
		total = c_add(total, *t->content->value);
	}

	return(total);
}


  
/****************************************************************************************************
*****************************************************************************************************
*****************************************************************************************************
*****************************************************************************************************
*****************************************************************************************************
*****************************************************************************************************/


