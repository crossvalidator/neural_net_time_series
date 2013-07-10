//#include "ancfis.h"
//#include "chaossimulatedannealling.c"

/*
This file is responsible for updating the parameters in layer 1.
It is also at this point where the constraints on the membership function
c & d are set.  The constraints so far are:

0<=c<=1
|d|<=c
|d|+c<=1

There are 2 functions:
1. update_parameter
	2. gradient_vector_length
*/

/***********************************************
Description:
Update parameters in k-th layer

Inputs:
k - the layer of the parameters to update. In this case
	k should be 1 since layer 4 parameters are only updated
	by the Kalman filter.  Also if k == -1, update all parameters
step_size - the step size

Outputs:
Updates the parameters using node_p[i]->parameter->content

Called from:
mexFunction in ancfismex.c
************************************************/
void update_parameter(int k, double step_size) {
	int i;
	int m=0;
	int j=0;
	int r = 0;

	double length, *magnitude, *phase, *chaos_SA;

	COMPLEX_LIST_T *p;
	COMPLEX_T tmp,*tmp1,tmp2;
	PARAMETER_LIST_T * plt;
	double gradient_vector_length();
	
	magnitude = calloc(In_vect_n*In_vect_n , sizeof(double));
	phase = calloc(In_vect_n*In_vect_n, sizeof(double));
	tmp1 = calloc(In_vect_n*In_vect_n , sizeof(COMPLEX_T));
   	length = gradient_vector_length(k);
	if (length == 0) {
		printf("gradient vector length == 0!\n");
		return;
	}
	
	for (i = 0; i < Node_n; i++) {
		if (node_p[i]->dE_dSMF == NULL)
			continue;
		if ((k != -1) && (node_p[i]->layer != k))
			continue;
		
    	tmp1 = first_initialized_weights(i);
    	
            
    	j = 0;
    	m = 0;
		for (p = node_p[i]->dE_dSMF; p != NULL; p = p->next) {	
		  	tmp = c_mul_scalar(*p->content, step_size/length);
		  	tmp2 = c_sub(tmp1[j++], tmp);
    	    magnitude[m] = c_abs(tmp2);
    	    phase[m] = c_phase(tmp2);
    	    m = m+1;
        }


		 chaos_SA = (double *)chaotic_simulatedAnnealling(magnitude, phase);
        r = 0;
        for (plt = node_p[i]->parameter ; plt != NULL ; plt = plt->next) {
            plt->content = chaos_SA[r++];
            //printf("parameter is :%f\n", plt->content);
        }
        
	}
	
}

/***********************************************
Description:
Calculate length of gradient vector of parameters in k-th layer.

Inputs:
k - the layer
	if k == -1, calculate length of gradient vector of all parameters

Outputs:
Returns the magnitude of the complex gradient vector length.

Called from:
update_parameter in new_para.c
************************************************/
double gradient_vector_length(int k) {
	int i;
	double total = 0;
	COMPLEX_LIST_T *p;
	for (i = 0; i < Node_n; i++) {
		if (node_p[i]->dE_dSMF == NULL)
			continue;
		if ((k != -1) && (node_p[i]->layer != k))
			continue;
		//moshkel dar node_p[i]->dE_dSMF ast ke adadhash kheiliajiband, chek konam 
		for (p = node_p[i]->dE_dSMF; p != NULL; p = p->next) {
		total += (p->content->real)*(p->content->real) + (p->content->imag)*(p->content->imag);
    	}
	}
	return(sqrt(total));
}
//*********************************************/

COMPLEX_T * first_initialized_weights(int node_index)
{
	PARAMETER_LIST_T *para_p = node_p[node_index]->parameter;
	double a, b, c, d, t, r, theta;
	int i, j, n=0;
	
	COMPLEX_T *weights;
	a = para_p->content;
	b = para_p->next->content;
	c = para_p->next->next->content;
	d = para_p->next->next->next->content;
	t = 2*PI/In_vect_n;
	
	weights = (COMPLEX_T *)create_array(In_vect_n*In_vect_n, sizeof(COMPLEX_T));
		
	for(i = 0 ; i < In_vect_n ; i++) {
	    theta = t*i;
	    r = d*sin(a*theta+b)+c;
	    
		for(j = 0 ; j < In_vect_n ; j++) { 
			weights[n++] = complex(r*cos(theta), r*sin(theta));		
		}
	}
   
	return(weights);	
}


