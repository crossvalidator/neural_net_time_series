//#include "ancfis.h"
/*
This file contains the code for the Kalman filter algorithm used
in layer 4.  The kalman filter is used in the forward pass (after calculating layer 3's output)
to determine the parameters {p, r}.

There are 11 functions:
1. get_kalman_data
2. kalman
	3. ma_plus_ma
	4. v_plus_v
	5. m_mult_v
	6. v_mult_m
	7. v_cross_v
	8. v_dot_v
	9. s_mult_ma
	10. s_mult_v
11. put_kalman_parameter
*/

static int k_p_n; //number of parameters to be identified by kalman filter algorithm

/***********************************************
Description:
Get the data to be used in kalman filter algorithm.
The kalman filter is an approximation of the least squares algorithm
kalman_data is the output

Inputs:
*kalman_data - array initialized in ancfismex.c
target - the last value in node_p

Outputs:
*kalman_data is filled with the data from node_p

			3 mf	layer 4 (DP) output
			/			a
x=[x1,x2,x3]-			b
			\			c

			/			d
y=[y1,y2,y3]-			e
			\			f

*kalman_data =  [ax1, ax2, ax3, ay1, ay2, ay3, a, bx1, bx2, bx3, by1, ... fy2, fy3, f, target]

Called from:
mexFunction in ancfismex.c
************************************************/
void get_kalman_data(double **kalman_data, double *target) {
	int h, i, k, m, j = 0;
	int first = In_n + In_n*Mf_n + 2*Rule_n; //'first' is the index of the first node in layer 4

	k_p_n = (In_vect_n*In_n +1)*Rule_n;


	for(m = 0 ; m < In_n ; m++) 
	{
		for(i = first ; i < first + Rule_n ; i++) 
		{
			for(k = 0 ; k < In_n ; k++) 
			{
				double *trn_data = node_p[k]->input_vector;
				for(h = 0 ; h < In_vect_n ; h++) 
				{
					kalman_data[m][j++] = (node_p[i]->value->real)*(trn_data[h]);
				}
			}
			kalman_data[m][j++] = node_p[i]->value->real;
		}
		j=0;
	}
}


/***********************************************
Description:
Sequential identification using kalman filter algorithm.
Note that s[][] and p[] are static variables, such that
their values are retained between every two function invocations.

Inputs:
zz - (zz+1) is the epoch number
k - (k+1) is the times 'kalman' has been invoked in an epoch
*kalman_data - from get_kalman_data
*kalman_parameter - the identified parameters

Outputs:
kalman_parameter stores the output values.

Called from:
mexFunction in ancfismex,c
************************************************/
void kalman(int zz, int k, int m, double **kalman_data, double **kalman_parameter, double *target) {
	void s_mult_v();
	void s_mult_ma();
	double v_dot_v();
	void v_cross_v();
	void v_mult_m();
	void m_mult_v();
	void v_plus_v();
	void ma_plus_ma();

	int i, j ;
	double alpha = 1000000.0;
	double denom, diff;
	static double **s;
	static double *p;
	static double *x, y;
	static double *tmp1, *tmp2;
	static double **tmp_m;

	//initial values of s[][] and p[] if k = 0

	if((zz == 0) && (k == 0)) {
		s = (double **)create_matrix(k_p_n, k_p_n, sizeof(double));
		p = (double *)create_array(k_p_n, sizeof(double));
		x = (double *)create_array(k_p_n, sizeof(double));
		tmp1 = (double *)create_array(k_p_n, sizeof(double));
		tmp2 = (double *)create_array(k_p_n, sizeof(double));
		tmp_m = (double **)create_matrix(k_p_n, k_p_n, sizeof(double));
	}

	if(k == 0) {
		for(i = 0; i < k_p_n; i++)
			p[i] = 0;

		for(i = 0; i < k_p_n; i++)
			for(j = 0; j < k_p_n; j++)
				if(i == j)
					s[i][j] = alpha;
				else
					s[i][j] = 0;
	}

	for(j = 0; j < k_p_n; j++)
		x[j] = kalman_data[m][j];
	y= target[m]; //target

	v_mult_m(x, s, tmp1);
	denom = 1 + v_dot_v(tmp1, x);
	m_mult_v(s, x, tmp1);
	v_mult_m(x, s, tmp2);
	v_cross_v(tmp1, tmp2, tmp_m);
	s_mult_ma(-1/denom, tmp_m);
	ma_plus_ma(s, tmp_m, s);

	diff = y - v_dot_v(x, p);
	m_mult_v(s, x, tmp1);
	s_mult_v(diff, tmp1);
	v_plus_v(p, tmp1, p);

	for(j = 0; j < k_p_n; j++) 
		kalman_parameter[m][j] = p[j];

}

/***********************************************
Description:
Matrix plus matrix.

Inputs:
**m1 - a matrix
**m2 - a matrix
**out - the result of m1 + m2

Outputs:
out is where the output is stored.

Called from:
kalman in kalman.c
************************************************/
void ma_plus_ma(double **m1, double **m2, double **out) {
	int i, j;
	for(i = 0; i < k_p_n; i++)
		for(j = 0; j < k_p_n; j++)
			out[i][j] = m1[i][j]+m2[i][j];
}

/***********************************************
Description:
Vector plus vector.

Inputs:
v1 - a vector
v2 - a vector
out - the result

Outputs:
out is where the output is stored.

Called from:
kalman in kalman.c
************************************************/
void v_plus_v(double v1[], double v2[], double out[]) {
	int i;
	for(i = 0; i < k_p_n; i++)
		out[i] = v1[i]+v2[i];
}

/***********************************************
Description:
Matrix multiplies vector.

Inputs:
**m - matrix
v[] - vector
out[] - output

Outputs:
out is where the output is stored.

Called from:
kalman in kalman.c
************************************************/
void m_mult_v(double **m, double v[], double out[]) {
	int i, j;
	for(i = 0; i < k_p_n; i++) {
		out[i] = 0;
		for(j = 0; j < k_p_n; j++)
			out[i] += m[i][j]*v[j];
	}
}

/***********************************************
Description:
Vector multiplies matrix

Inputs:
v[] - vector
**m - matrix
out[] - output

Outputs:
out is where the output is stored.

Called from:
kalman in kalman.c
************************************************/
void v_mult_m(double v[], double **m, double out[]) {
	int i, j;
	for(i = 0; i < k_p_n; i++) {
		out[i] = 0;
		for(j = 0; j < k_p_n; j++)
			out[i] += v[j]*m[j][i];
	}
}


/***********************************************
Description:
Vector cross-product

Inputs:
v1[] - vector
v2[] - vector
**out - output

Outputs:
out is where the output is stored.

Called from:
kalman in kalman.c
************************************************/
void v_cross_v(double v1[], double v2[], double **out) {
	int i, j;
	for(i = 0; i < k_p_n; i++)
		for(j = 0; j < k_p_n; j++)
			out[i][j] = v1[i]*v2[j];
}

/***********************************************
Description:
Vector dot product

Inputs:
v1[] - vector
v2[] - vector

Outputs:
Returns the dot product.

Called from:
kalman in kalman.c
************************************************/
double v_dot_v(double v1[], double v2[]) {
	int i;
	double out = 0;
	for(i = 0; i < k_p_n; i++)
		out += v1[i]*v2[i];
	return(out);
}

/***********************************************
Description:
Scalar multiplies matrix

Inputs:
c - scalar value
**m - matrix

Outputs:
Changes stay in matrix m.

Called from:

************************************************/
void s_mult_ma(double c, double **m) {
	int i, j;
	for(i = 0; i < k_p_n; i++)
		for(j = 0; j < k_p_n; j++)
			m[i][j] = c*m[i][j];
}

/***********************************************
Description:
Scalar multiplies vector

Inputs:
c - scalar value
v[] - vector

Outputs:
Changes stay in vector v.

Called from:

************************************************/
void s_mult_v(double c, double v[]) {
	int i;
	for(i = 0; i < k_p_n; i++)
		v[i] = c*v[i];
}

/***********************************************
Description:
Put the parameters identified by kalman filter algorithm into GNN.
The parameters in our case are {p, r}.

Inputs:
*kalman_parameter - the identified parameters

Outputs:
Operates on node_p[i]->parameter for layer 4

Called from:
mexFunction in ancfismex.c
************************************************/
void put_kalman_parameter(int Out_n,double **kalman_parameter) {
	//'first' is the index of the first node in layer 5
	int first = In_n + In_n*Mf_n + 3*Rule_n;
	int i,k, j;
	PARAMETER_LIST_T *p;

	for (k=0; k < Out_n; k++)
	{
		j=0;
		for(i = first+k; i < first + (Rule_n * Out_n) ; i=i+Out_n)
		{
			for(p = node_p[i]->parameter; p != NULL; p = p->next)
				p->content = kalman_parameter[k][j++];
		}
	}
}


