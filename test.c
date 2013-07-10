#include "mex.h"
#include "ancfis.h"
#include "lib.c"
#include "complex.c"

void mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
const mxArray *prhs[];
mxArray *plhs[];
{
    int i = 0, j = 0;
	int In_vect_n = 3;    
    COMPLEX_T* a;
    COMPLEX_T tmp;
   
    for(j = 0 ; j < In_vect_n ; j++) {
        a = first_initialized_weights(j);
    
        for(i = 0 ; i < In_vect_n*In_vect_n ; i++) {
            tmp = a[i];
            mexPrintf("%f, %f\n", tmp.real, tmp.imag);
        }
    }
}

COMPLEX_T* first_initialized_weights(int k)
{
	double a, b, c, d, t, r, theta;
	int i, j, n = 0;
	int In_vect_n = 3;
	
	COMPLEX_T *weights;
	a = 0.1*k;
	b = 0.2*k;
	c = 0.3*k;
	d = 0.4*k;
	t = 2*PI/In_vect_n;
	
	weights = (COMPLEX_T *)create_array(In_vect_n*In_vect_n, sizeof(COMPLEX_T));
		
	for(i = 0 ; i < In_vect_n ; i++) {
	    theta = t*i;
	    r = d*sin(a*theta+b)+c;
	    
		for(j = 0 ; j < In_vect_n ; j++) { 
			weights[n++] = complex(r*cos(theta), r*sin(theta));
			mexPrintf("weight=%f ,%f \n", weights[n-1].real, weights[n-1].imag);		
		}
	}
   
	return(weights);
}
