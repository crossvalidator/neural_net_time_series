
// please read only this instruction.
// this is new code after discussion with prof scott dick. here, layer 3 and alyer 4 are redefined.
// this is second improveemnt
//#include "ancfis.h"
// this is a basic code of teh original entire code... many other improvements have been made in this code and these 
// improved backward.c codes are typed below this code with proper commenting. Please have a look at them too.
/*
This file calculates the error signals (derivatives) in the backward pass.
This excludes the derivatives with respect to each parameter
as those are done in the file de_dp.

There are 7 functions:
1. calculate_de_do
	2. derivative_o_o
	    3. do6_do5
		4. do5_do4
		5. do4_do3
		6. do3_do2
		7. do2_do1
*/

/***********************************************
Description:
This calculates de/do of each node. de/do is the (e)rror wrt the (o)utput

Note that de/do = dE/dx = epsilon in Zhifei's paper

Inputs:
de_dout - is de/do of the very last node in the network.
	For sum of squares error, de_dout = -2*(actual data - predicted output).

Outputs:
N/A

Called from:
mexFunction in ancfismex.c
*********************************************************************************/

void calculate_de_do(int k,double Out_n,double *target,double **ancfis_output) {
	int i, j;
		//,m;
	NODE_LIST_T *p;
	COMPLEX_T tmp1, tmp2;
	COMPLEX_T derivative_o_o();
	COMPLEX_T de_do;
	COMPLEX_T de_dout;

	for (i=0; i<Out_n; i++)
	{
		//get the desired output
		de_dout = complex(target[i] - ancfis_output[k][i], 0); 
		//Calculate the error rate for output node
		de_dout = c_mul_scalar(de_dout, -2.0);
	//	m= (int)(Out_n-i);
		*node_p[Node_n - (i+1)]->de_do = de_dout; //set the very last node
	}

	for(i = Node_n - (int)(Out_n+1); i >= In_n; i--) { //start from the 2nd last node
		de_do = complex(0.0, 0.0);

			for(p = node_p[i]->fan_out; p != NULL; p = p->next) {
			j = p->content->index; //index of some node in the next layer
			tmp1 = *p->content->de_do; //de_do of some node in the next layer
			tmp2 = derivative_o_o(i, j);
			tmp1 = c_mul(tmp1, tmp2);
			de_do = c_add(de_do, tmp1);
		}

		*node_p[i]->de_do = de_do;
	}
}

/***********************************************
Description:
Calculate do(i)/do(j), where i and j are node indices.
This function simply calls other functions based on the layer.
do/do is the derivative of an (o)utput wrt the (o)utput of another layer.

Note: do/do is equvalent to df/dx in Zhifei's paper

Inputs:
i and j - node indices

Outputs:
Returns do(i)/do(j).

Called from:
calculate_de_do in backward.c
************************************************/

COMPLEX_T derivative_o_o(int i, int j) {
	COMPLEX_T do2_do1(), do3_do2(), do4_do3(), do5_do4(),do6_do5() ;
	int layer = node_p[i]->layer;

	switch(layer) {
		case 0:
			printf("Error in derivative_o_o!");
		case 1:
			return(do2_do1(i, j));
		case 2:
			return(do3_do2(i, j));
		case 3:
			return(do4_do3(i, j));
		case 4:
			return(do5_do4(i, j));
	    case 5:
	        return(do6_do5(i, j));
		default:
			printf("Error in derivative_o_o!");
	}

	return(complex(0.0, 0.0)); //suppress compiler error
}

/***********************************************
Description:
Calculate do(i)/do(j), where node i is in layer 5, node j is in layer 6
For our network, the derivative from layer 6 to 5 is 1.

                        0
                        |\
     mf1--3-----5-----7---9
   /       \   / \   /  |   \
 0           \     /    |     11--
   \       /   \ /   \   \  /
     mf2--4-----6-----8---10

input 1   2     3     4   5 (layer)

Inputs:
i - node index in layer 5
j - node index in layer 6

Outputs:
Returns do(i)/do(j) = 1.0.

Called from:
derivative_o_o in backward.c
************************************************/

COMPLEX_T do6_do5(int i, int j) {
	return(complex(1.0, 0.0));
}

/***********************************************
Description:
Calculate do(i)/do(j), where node i is in layer 4, node j is in layer 5.

From the forward pass, the output of node 7 is W1dp
The output of node 9 is W1dp*(p1x1+p2x2+ ... + r)
The resulting derivative is simply (p1x1+p2x2+ ... + r)

                        0
                        |\
     mf1--3-----5-----7---9
   /       \   / \   /  |   \
 0           \     /    |     11--
   \       /   \ /   \   \  /
     mf2--4-----6-----8---10

input 1   2     3     4   5 (layer)

Inputs:
i - node index in layer 4
j - node index in layer 5

Outputs:
Returns do(i)/do(j) = (p1x1+p2x2+ ... + r)

Called from:
derivative_o_o in backward.c
************************************************/

COMPLEX_T do5_do4(int i, int j) {
	PARAMETER_LIST_T *para_p = node_p[j]->parameter;
	PARAMETER_LIST_T *p_p;
	int k, h;
	double a, tmp, total = 0.0;

	p_p = para_p;

	for(k = 0; k < In_n + 1; k++) {
		double *trn_data = node_p[k]->input_vector;
		if(k == In_n)
			break;

		for(h = 0 ; h < In_vect_n ; h++) {
            a = p_p->content;//this is p1, p2 ... q1, q2...
			tmp = a*(trn_data[h]); //p1x1, p2x2, ...
			total += tmp;
			p_p = p_p->next;
		}
	}

	a = p_p->content; //this is r
    tmp = total+a;
	return(complex(tmp, 0.0));

}

/***********************************************
Description:
Calculate do(i)/do(j), where node i is in layer 3, node j layer 4.

                        0
                        |\
     mf1--3-----5-----7---9
   /       \   / \   /  |   \
 0           \     /    |     11--
   \       /   \ /   \   \  /
     mf2--4-----6-----8---10

input 1   2     3     4   5 (layer)

Inputs:
i - node index in layer 3
j - node index in layer 4

Outputs:
Returns do(i)/do(j).

Called from:
derivative_o_o in backward.c
************************************************/
/*
COMPLEX_T do4_do3(int i, int j) {
	NODE_LIST_T *arg_p = node_p[j]->fan_in;
	NODE_LIST_T *p;
	//COMPLEX_T total = complex(0.0, 0.0);
	COMPLEX_T tmp;
//, tmp2, conjg_total;
/*
	for(p = arg_p; p != NULL; p = p->next) {
		total = c_add(total, *p->content->value); //sum of previous layer
	}
*/
	//conjg_total = c_conjg(total);
/*
	if((j - i) == Rule_n) { //Rule_n is # of nodes in 4th layer
		//straight connection
		//tmp = c_mul_scalar(conjg_total, 0.5);
		//tmp2 = c_conjg(*node_p[i]->value);
		//tmp2 = c_mul_scalar(tmp2, 0.5);
		//tmp2 = c_add(tmp, tmp2);
		tmp = complex(0.5, 0.0);
		return(tmp);
	} else {
		//cross connection
		tmp = complex(-0.5, 0.0);
		return(tmp);
	}
}
*/
/************************************************
COMPLEX_T do4_do3(int i, int j) {
	NODE_LIST_T *arg_p = node_p[j]->fan_in;
	NODE_LIST_T *p;
	COMPLEX_T total = complex(0.0, 0.0);
	COMPLEX_T tmp, tmp2, conjg_total;

	for(p = arg_p; p != NULL; p = p->next) {
		total = c_add(total, *p->content->value); //sum of previous layer
	}

	conjg_total = c_conjg(total);

	if((j - i) == Rule_n) { //Rule_n is # of nodes in 4th layer
		//straight connection
		tmp = c_mul_scalar(conjg_total, 0.5);
		tmp2 = c_conjg(*node_p[i]->value);
		tmp2 = c_mul_scalar(tmp2, 0.5);
		tmp2 = c_add(tmp, tmp2);
		return(tmp2);
	} else {
		//cross connection
		tmp = c_conjg(*node_p[j - Rule_n]->value);
		return(c_mul_scalar(tmp, 0.5));
	}
}
***************************************/

COMPLEX_T do4_do3(int i, int j) {
	NODE_LIST_T *arg_p = node_p[j]->fan_in;
	NODE_LIST_T *p;
	//COMPLEX_T total = complex(0.0, 0.0);
	COMPLEX_T temp1, temp2;
	COMPLEX_T gen_sum = complex(0.0, 0.0);
	//double abs_gen_sum = 0.0;
	//double sqr_mag_sum = 0.0;
	//double sqrof_abs_gen_sum = 0.0;
	//double complicated_diff = 0.0;
	//double sqrt_complicated_diff = 0.0;

	for(p = arg_p; p != NULL; p = p->next) { //loop for all fan_in nodes
		//mag_sum += c_abs(*p->content->value);
		gen_sum = c_add(gen_sum, *p->content->value);
              	//sqr_mag_sum = sqr_mag_sum + (c_abs(*p->content->value)*c_abs(*p->content->value));
	}

	//abs_gen_sum = c_abs(gen_sum);
	//sqrof_abs_gen_sum = abs_gen_sum * abs_gen_sum;
	//complicated_diff = sqrof_abs_gen_sum - sqr_mag_sum;
	//sqrt_complicated_diff = sqrt(fabs(complicated_diff));


if((j - i) == Rule_n) { //Rule_n is # of nodes in 4th layer
		//straight connection

//temp1= c_conjg(gen_sum);
temp2= *node_p[i]->value;
temp1= c_conjg(c_sub(gen_sum, temp2));
return(c_mul_scalar(temp1, 1.0/2.0));
}

else{
temp2= *node_p[j-Rule_n]->value;
temp1= c_conjg(temp2);
return(c_mul_scalar(temp1, 1.0/2.0));
}
}




//int test_total = c_abs(total);



/*********************************************
if(test_total !=0){
	if((j - i) == Rule_n) { //Rule_n is # of nodes in 4th layer
		//straight connection
		tmp = c_mul_scalar(conjg_total, 0.5);
		tmp2 = c_conjg(*node_p[i]->value);
		tmp2 = c_mul_scalar(tmp2, 0.5);
		tmp2 = c_add(tmp, tmp2);
		return(tmp2);
			       } 
				else {
					//cross connection
					tmp = c_conjg(*node_p[j - Rule_n]->value);
					return(c_mul_scalar(tmp, 0.5));
				     }
		  }

	else {
		if((j - i) == Rule_n) 
				{ 
					//Rule_n is # of nodes in 4th layer
					//straight connection
					//tmp = c_mul_scalar(conjg_total, 0.5);
					tmp2 = c_conjg(*node_p[i]->value);
					//tmp2 = c_mul_scalar(tmp2, 0.5);
					//tmp2 = c_add(tmp, tmp2);
					return(tmp2);
				}
			else
			{
				return (complex(0.0, 0.0));
			}

	     }



***********************************/


	


/***********************************************
Description:
Calculate do(i)/do(j), where node i is in layer 2, node j layer 3.

                        0
                        |\
     mf1--3-----5-----7---9
   /       \   / \   /  |   \
 0           \     /    |     11--
   \       /   \ /   \   \  /
     mf2--4-----6-----8---10

input 1   2     3     4   5 (layer)

Inputs:
i - node index in layer 2
j - node index in layer 3

Outputs:
Returns do(i)/do(j).

Called from:
derivative_o_o in backward.c
************************************************/

COMPLEX_T do3_do2(int i, int j) {
	NODE_LIST_T *arg_p = node_p[j]->fan_in;
	NODE_LIST_T *p;
	double ux, uy, vx, vy, mag_sum = 0.0;
	COMPLEX_T actual_sum;
	actual_sum = complex(0.0, 0.0);
	COMPLEX_T temp, temp2;
	for(p = arg_p; p != NULL; p = p->next) { //loop for all fan_in nodes
		actual_sum = c_add(actual_sum, *p->content->value);
	}
	
	if((j - i) == Rule_n) 
	{
	temp = *node_p[i]->value;
	temp = c_sub(actual_sum, temp);
	temp2 = c_mul(actual_sum,actual_sum);
	temp = c_div(temp, temp2);
	return(temp);
	}
	else
	{
	temp = *node_p[j-Rule_n]->value;
	temp = c_mul_scalar(temp, -1.0);
	temp2 = c_mul(actual_sum,actual_sum);
	temp = c_div(temp, temp2);
	return(temp);
	}
}

/******************************************
	if((j - i) == Rule_n) { //Rule_n is # of nodes in 4th layer
		//straight connection
		ux = (mag_sum - ((node_p[i]->value->real*node_p[i]->value->real)/c_abs(*node_p[i]->value)))/(mag_sum*mag_sum);
		uy = (-node_p[i]->value->real*node_p[i]->value->imag)/c_abs(*node_p[i]->value)/(mag_sum*mag_sum);
		vx = uy;
		vy = (mag_sum - ((node_p[i]->value->imag*node_p[i]->value->imag)/c_abs(*node_p[i]->value)))/(mag_sum*mag_sum);
	} else {
		//cross connection
		ux = (-node_p[j - Rule_n]->value->real*node_p[i]->value->real)/c_abs(*node_p[i]->value)/(mag_sum*mag_sum);
		uy = (-node_p[j - Rule_n]->value->real*node_p[i]->value->imag)/c_abs(*node_p[i]->value)/(mag_sum*mag_sum);
		vx = (-node_p[j - Rule_n]->value->imag*node_p[i]->value->real)/c_abs(*node_p[i]->value)/(mag_sum*mag_sum);
		vy = (-node_p[j - Rule_n]->value->imag*node_p[i]->value->imag)/c_abs(*node_p[i]->value)/(mag_sum*mag_sum);
	}

	return(complex(0.5*(ux+vy), 0.5*(-uy+vx)));
**********************************/



/***********************************************
Description:
Calculate do(i)/do(j), where node i is in layer 1, node j layer 2.

                        0
                        |\
     mf1--3-----5-----7---9
   /       \   / \   /  |   \
 0           \     /    |     11--
   \       /   \ /   \   \  /
     mf2--4-----6-----8---10

input 1   2     3     4   5 (layer)

Inputs:
i - node index in layer 1
j - node index in layer 2

Outputs:
Returns do(i)/do(j).

Called from:
derivative_o_o in backward.c
************************************************/

COMPLEX_T do2_do1(int i, int j) {
	return(c_div(*node_p[j]->value, *node_p[i]->value));
}


/**************************************************************************************************
***************************************************************************************************
***************************************************************************************************
***************************************************************************************************
***************************************************************************************************
**************************************************************************************************/


