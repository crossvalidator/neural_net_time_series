//#include "ancfis.h"
#define PI 3.1415926535897932385E0
/*
These are functions for performing complex arithmetic.

There are 12 functions:
1. complex
2. c_add
3. c_sub
4. c_mul
5. c_conjg
6. c_div
7. c_abs
8. c_sqrt
9. c_mul_scalar
10. c_add_scalar
11. c_phase
12. c_dot_product
*/

/***********************************************
Description:
This returns a new COMPLEX_T number.  This is
like a constructor.

Inputs:
re - the real portion
im - the imaginary portion

Outputs:
A new COMPLEX_T number.
************************************************/
COMPLEX_T complex(double re, double im) {
	COMPLEX_T c;
	c.real = re;
	c.imag = im;
	return c;
}

/***********************************************
Description:
This is complex addition a+b.

Inputs:
a - the first complex number to be added
b - the second complex number to be added

Outputs:
a+b
************************************************/
COMPLEX_T c_add(COMPLEX_T a, COMPLEX_T b) {
	COMPLEX_T c;
	c.real = a.real+b.real;
	c.imag = a.imag+b.imag;
	return c;
}

/***********************************************
Description:
This is complex subtraction a-b.

Inputs:
a - the first complex number
b - the complex number that is subtracted from the first

Outputs:
a-b
************************************************/
COMPLEX_T c_sub(COMPLEX_T a, COMPLEX_T b) {
	COMPLEX_T c;
	c.real = a.real-b.real;
	c.imag = a.imag-b.imag;
	return c;
}

/***********************************************
Description:
This is complex multiplication a*b.

Inputs:
a - the first complex number to be multiplied
b - the second complex number to be multiplied

Outputs:
a*b
************************************************/
COMPLEX_T c_mul(COMPLEX_T a, COMPLEX_T b) {
	COMPLEX_T c;
	if (fabs(a.real) >= (1e+010)) (a.real=0);
	if (fabs(b.real) >= (1e+010)) (b.real=0);
	if (fabs(a.imag) >= (1e+010)) (a.imag=0);
	if (fabs(b.imag) >= (1e+010)) (b.imag=0);
	c.real = a.real*b.real-a.imag*b.imag;
	c.imag = a.imag*b.real+a.real*b.imag;
	return c;
}

/***********************************************
Description:
Returns the complex conjugate of z = x + yj which is
x - yj.

Inputs:
z - a complex number

Outputs:
The congugate of z
************************************************/
COMPLEX_T c_conjg(COMPLEX_T z) {
	COMPLEX_T c;
	c.real = z.real;
	c.imag = -z.imag;
	return c;
}

/***********************************************
Description:
Returns the division of complex numbers a/b.

Inputs:
a - numerator
b - denominator

Outputs:
a/b
************************************************/
COMPLEX_T c_div(COMPLEX_T a, COMPLEX_T b) {
	COMPLEX_T c;
	double r,den;
	if (fabs(b.real) >= fabs(b.imag)) {
		r = b.imag/b.real;
		den = b.real+r*b.imag;
		c.real = (a.real+r*a.imag)/den;
		c.imag = (a.imag-r*a.real)/den;
	} else {
		r = b.real/b.imag;
		den = b.imag+r*b.real;
		c.real = (a.real*r+a.imag)/den;
		c.imag = (a.imag*r-a.real)/den;
	}
	return c;
}

/***********************************************
Description:
Returns the absolute value or the magnitude of a complex number.

Inputs:
z - a complex number

Outputs:
sqrt(x^2 + y^2)
************************************************/
double c_abs(COMPLEX_T z) {
	/*double x,y,ans,temp;
	x = fabs(z.real);
	y = fabs(z.imag);
	if (x == 0.0)
		ans = y;
	else if (y == 0.0)
		ans = x;
	else if (x > y) {
		temp = y/x;
		ans = x*sqrt(1.0+temp*temp);
	} else {
		temp = x/y;
		ans = y*sqrt(1.0+temp*temp);
	}
	return ans;*/

	//mexPrintf("'%f, %f'\n",z.real,z.imag);
	return sqrt(z.real*z.real + z.imag*z.imag);
}

/***********************************************
Description:
Returns the square root of a complex number.
+	z	{...}

Inputs:
z - a complex number

Outputs:
sqrt(z)
************************************************/
COMPLEX_T c_sqrt(COMPLEX_T z) {
	COMPLEX_T c;
	double x,y,w,r;
	if ((z.real == 0.0) && (z.imag == 0.0)) {
		c.real = 0.0;
		c.imag = 0.0;
		return c;
	} else {
		x = fabs(z.real);
		y = fabs(z.imag);
		if (x >= y) {
			r = y/x;
			w = sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r = x/y;
			w = sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.real >= 0.0) {
			c.real = w;
			c.imag = z.imag/(2.0*w);
		} else {
			c.imag = (z.imag >= 0) ? w : -w;
			c.real = z.imag/(2.0*c.imag);
		}
		return c;
	}
}

/***********************************************
Description:
Multiply a complex number by a scalar.  Both real and imaginary
parts are scaled.

Inputs:
a - a complex number
x - a scalar

Outputs:
a*x
************************************************/
COMPLEX_T c_mul_scalar(COMPLEX_T a, double x) {
	COMPLEX_T c;
	c.real = x*a.real;
	c.imag = x*a.imag;
	return c;
}

/***********************************************
Description:
Add a complex number to a scalar.  Only the real portion is affected.

Inputs:
a - a complex number
x - a scalar

Outputs:
a+x
************************************************/
COMPLEX_T c_add_scalar(COMPLEX_T a, double x) {
	COMPLEX_T c;
	c.real = x+a.real;
	c.imag = a.imag;
	return c;
}

/***********************************************
Description:
Returns the phase of a complex number.

Inputs:
z - a complex number

Outputs:
atan(fabs(y)/fabs(x)) or PI-atan(fabs(y)/fabs(x)) or PI+atan(fabs(y)/fabs(x)) or 2*PI-atan(fabs(y)/fabs(x))
************************************************/

double c_phase(COMPLEX_T z) {
	double x,y,ans;
	x = z.real;
	y = z.imag;
	if (x >= 0.0){
		if (y >= 0)
		    ans = atan(fabs(y)/fabs(x));
		else
		    ans = 2*PI-atan(fabs(y)/fabs(x));}
	else {
		if (y >= 0.0)
		    ans = PI-atan(fabs(y)/fabs(x));
		else
		    ans = PI+atan(fabs(y)/fabs(x));}
	return ans;
}

/***********************************************
Description:
This returns the dot product a dot b.  Note that the complex number is treated as a
regular vector so no conjugate is needed.

It returns a complex number with the imaginary part set to 0.

Inputs:
a - the first complex number
b - the second complex number

Outputs:
a dot b
************************************************/
COMPLEX_T c_dot_product(COMPLEX_T a, COMPLEX_T b) {
	COMPLEX_T c;
	c.real = a.real*b.real + a.imag*b.imag;
	c.imag = 0;
	return c;
}
