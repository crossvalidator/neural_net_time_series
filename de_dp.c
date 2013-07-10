//#include "ancfis.h"
/*
This file contains functions that calculate the error derivatives
with respect to the parameters.

In our case, only layers 1 and 4 have parameters.

There are 5 functions:
1. clear_de_dp
2. update_de_dp
	3. derivative_o_p
		4. dmf_dp
		5. dconsequent_dp
*/

/***********************************************
Description:
Clear de/dp of each parameter.  Iterate through all of the nodes
and set de_dp = 0.  This needs to be done in between epochs since the
same data structure is being used for all epochs.

Inputs:
N/A

Outputs:
Operates on global variable 'node_p'

Called by:
mexFunction in ancfismex.c
************************************************/
void clear_de_dp() {
	int i;
	COMPLEX_LIST_T *p;
	COMPLEX_T tmp;

	for(i = 0; i < Node_n; i++) {
		if(node_p[i]->dE_dSMF == NULL)
			continue;
		for(p = node_p[i]->dE_dSMF; p != NULL; p = p->next) {
			tmp = complex(0.0, 0.0);
			*p->content = tmp;
		}
	}
}

/***********************************************
Description:
Calculate de/dp of each parameter.

In Zhifei's paper, de/dp = dE/SMF

Inputs:
iteration_number - the number of the current training data pair

Outputs:
Operates on global variable 'node_p'

Called by:
mexFunction in ancfismex.c
************************************************/
void update_de_dp(int iteration_number) {
	int i, j;
	COMPLEX_LIST_T *p;
	COMPLEX_T do_dp, tmp, tmp2;
	COMPLEX_T derivative_o_p();

	for(i = 0; i < Node_n; i++) {
		if(node_p[i]->dE_dSMF == NULL)
			continue;

		j = 0;

		for(p = node_p[i]->dE_dSMF; p != NULL; p = p->next) {
			do_dp = derivative_o_p(i, j, iteration_number);
			tmp = c_mul(*node_p[i]->de_do, do_dp);
			tmp2 = c_add(*p->content, tmp);
			*p->content = tmp2;
			j++;

			if(j == In_vect_n) {
			    j = 0;
			}
		}
	}
}

/***********************************************
Description:
Calculate the derivative of node i wrt it's j-th parameter.
Since only layer 1 and 5 have parameters, this function calls
other functions based on the layer.  Also since layer 5 is not
updated using this method, it is not called.

In Zhifei's paper, do/dp = df/SMF

Inputs:
i - node index
j - j-th parameter
iteration_number - the number of the current training data pair

Outputs:
The derivative of node i wrt it's j-th parameter.

Called by:
update_de_dp in de_dp.c
************************************************/
COMPLEX_T derivative_o_p(int i, int j, int iteration_number) {
	COMPLEX_T dmf_dp();
	COMPLEX_T dconsequent_dp();
	int layer = node_p[i]->layer;

	switch(layer) {
		case 1:
			return(dmf_dp(i, j, iteration_number));
		default:
			printf("Error in derivative_o_p!");
	}

	return(complex(0.0, 0.0));//suppress compiler warning
}

/***********************************************
Description:
Calculate the derivative of node i wrt it's j-th parameter
for layer 1.

Here, we have parameters {a,b,c,d} for the membership function:
-dsin(-(at+b))+c

Inputs:
i - node index
j - j-th parameter
iteration_number - the number of the current training data pair

Outputs:
The derivative of layer 1
-dsin(-(at+b))+c w.r.t. {a,b,c,d}

Called by:
derivative_o_p in de_dp.c
************************************************/
COMPLEX_T dmf_dp(int i, int j, int iteration_number) {
	NODE_LIST_T *arg_p = node_p[i]->fan_in ;
	double *input_vect = arg_p->content->input_vector;
	double ux, uy, vx, vy, abs_sum_conv, denom;
	COMPLEX_T sum_conv, deriv_sum_conv, df_dz;

	sum_conv = *node_p[i]->value;
	deriv_sum_conv = complex(input_vect[j], 0);
	abs_sum_conv = c_abs(sum_conv);

	if(abs_sum_conv == 0.0) {
		return deriv_sum_conv;
	} else {
		denom = abs_sum_conv*pow(1 + abs_sum_conv, 2);
		ux = (pow(sum_conv.imag, 2) + abs_sum_conv)/denom;
		uy = sum_conv.real*sum_conv.imag/(-denom);
		vx = uy;
		vy = (pow(sum_conv.real, 2) + abs_sum_conv)/denom;
	}

	df_dz = complex(0.5*(ux + vy), 0.5*(-uy + vx));
	return(c_mul(deriv_sum_conv, df_dz));
}


/***********************************************
Description:
Calculate the derivative of node i wrt it's j-th parameter
for layer 5.  Here we have {p,r} in w(px+r).  Note that although
this code exists, this is never actually used because the parameters
in this layer are updated using a LSE in the forward pass.

Inputs:
i - node index
j - node index

Outputs:
The derivative of layer 4 wrt its parameters {p,r}

Called by:
derivative_o_p in de_dp.c
************************************************/
COMPLEX_T dconsequent_dp(int i, int j) {
	NODE_LIST_T *arg_p = node_p[i]->fan_in ;
	PARAMETER_LIST_T *para_p = node_p[i]->parameter;
	NODE_LIST_T *a_p;
	int k;
	COMPLEX_T wn, x;

	//advance to the last fan in node
	for(a_p = arg_p; a_p->next != NULL; a_p = a_p->next);

	wn = *a_p->content->value;

	if(j == In_n)
		return(wn);

	for(a_p = arg_p, k = 0; k < j; a_p = a_p->next, k++);

	x = *a_p->content->value;

	return(c_mul(x, wn));
}
