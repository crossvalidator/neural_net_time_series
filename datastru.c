//#include "ancfis.h"
/*
This file contains functions that allocate memory for the
ANCFIS data structure.
There are 6 functions.
1. build_ancfis
	2. build_layer
		3. build_parameter_list
		4. build_complex_list
	5. build_node_list
6. set_parameter_mode
*/
/***********************************************
Description:
This function allocates memory/sets parameters for the ANCFIS data
structure by calling other functions.

Inputs:
**config - a square matrix generated by gen_config in gen_config.c.
	It is a matrix of dimensions (total number of nodes)x(total number of nodes).
	It contains only 0's and 1's which indicate if there is a forward connection
	between the nodes.  For example,

0 1
0 0

would show that node 0 is connected to node 1.

Outputs:
N/A

Called from:
mexFunction in ancfismex.c
************************************************/
void build_ancfis(int **config) {
	void build_layer();
	NODE_LIST_T *build_node_list();
	int i;
	build_layer(0, In_n, 0, 0, 0, 0);
	build_layer(1, In_n*Mf_n, In_n, 4, 1, In_vect_n*In_vect_n);// 4 is a,b,c,d
	build_layer(2, Rule_n, In_n+In_n*Mf_n, 0, 2, 0);
	build_layer(3, Rule_n, In_n+In_n*Mf_n+Rule_n, 0, 3, 0);
	build_layer(4, Rule_n, In_n+In_n*Mf_n+2*Rule_n,0, 4, 0);
	build_layer(5, Rule_n*Out_n, In_n+In_n*Mf_n+3*Rule_n, In_n*In_vect_n + 1, 5, 0);
	build_layer(6, Out_n, In_n+In_n*Mf_n+3*Rule_n+Rule_n*Out_n, 0, 6, 0);

	for(i = 0; i < Node_n; i++) {
		node_p[i]->fan_in  = build_node_list(0, i, config);
		node_p[i]->fan_out = build_node_list(1, i, config);
	}
}

/***********************************************
Description:
Build a node list of layer-th layer, with n nodes, starting at index index.
Each node has parameter_n parameters and node function
function[function_index]. Also calls a function to allocate memory for the
parameter list of the node.

Inputs:
layer - the layer (0-5)
n - number of nodes
index - starting node index
parameter_n - number of parameters in each node (4 for layer 1, In+1 for layer 4)
function_index - defines what function is used in the forward pass, (0-5)
dE_dSMF_n - number of complex values in dE_dSMF_n

Outputs:
Allocates memory/sets parameters for the global variable 'node_p'

Called from:
build_ancfis in datastru.c
************************************************/
void build_layer(int layer, int n, int index, int parameter_n, int function_index, int dE_dSMF_n) {
	PARAMETER_LIST_T *build_parameter_list();
	int i;
	NODE_T *q;

	if(n == 0)
		printf("Possible error in build_layer!\n");

	for(i = 0; i < n; i++)
	{
		q = malloc(sizeof(NODE_T));
		node_p[i + index] = q;
		q->index = i + index;
		q->layer = layer;
		q->local_index = i;
		q->parameter_n = parameter_n;
		q->function_index = function_index;

		q->value = malloc(sizeof(COMPLEX_T));
		q->de_do = malloc(sizeof(COMPLEX_T));
		q->value->real = 0;
		q->value->imag = 0;
		q->de_do->real = 0;
		q->de_do->imag = 0;

		q->dE_dSMF = build_complex_list(dE_dSMF_n);

		q->parameter = build_parameter_list(parameter_n);
	}
}

/***********************************************
Description:
Build a parameter list with length n.

Inputs:
n - the number of parameters

Outputs:
Returns the PARAMETER_LIST_T structure

Called by:
build_layer in datastru.c
************************************************/
PARAMETER_LIST_T * build_parameter_list(int n) {
	PARAMETER_LIST_T *head, *p, *q;
	int i;

	if(n < 0)
		printf("Error in build_parameter_list()!");

	if(n == 0)
		return(NULL);

	head = malloc(sizeof (PARAMETER_LIST_T));
	head->de_dp = malloc(sizeof(COMPLEX_T));
	head->de_dp->real = 0.0;
	head->de_dp->imag = 0.0;
	p = head;

	for(i = 1; i < n; i++){
		q = malloc(sizeof (PARAMETER_LIST_T));
		q->de_dp = malloc(sizeof(COMPLEX_T));
		q->de_dp->real = 0.0;
		q->de_dp->imag = 0.0;
		p->next = q;
		p = q;
	}

	p->next = NULL;
	return(head);
}

/***********************************************
Description:
Build a complex list with length n.

Inputs:
n - the number of complex values in the list

Outputs:
Returns the COMPLEX_LIST_T structure

Called by:
build_layer in datastru.c
************************************************/
COMPLEX_LIST_T * build_complex_list(int n) {
	COMPLEX_LIST_T *head, *p, *q;
	int i;

	if(n < 0)
		printf("Error in build_complex_list()!");

	if(n == 0)
		return(NULL);

	head = malloc(sizeof (COMPLEX_LIST_T));
	head->content = malloc(sizeof(COMPLEX_T));
	head->content->real = 0.0;
	head->content->imag = 0.0;
	p = head;

	for(i = 1; i < n; i++){
		q = malloc(sizeof (COMPLEX_LIST_T));
		q->content = malloc(sizeof(COMPLEX_T));
		q->content->real = 0.0;
		q->content->imag = 0.0;
		p->next = q;
		p = q;
	}

	p->next = NULL;
	return(head);
}

/***********************************************
Description:
Build the fan in and fan out node lists for each node.
Basically define the connections between the nodes.

Inputs:
type - (0,1)
	0 --> build node list along column n of matrix config.
	1 --> build node list along row n of matrix config.
n - node index
**config - a square matrix generated by gen_config in gen_config.c.
	It is a matrix of dimensions (total number of nodes)x(total number of nodes).
	It contains only 0's and 1's which indicate if there is a forward connection
	between the nodes.

Outputs:
NODE_LIST_T structure which is the fan in node list (type == 0)
or fan out node list (type == 1).

Called by:
build_ancfis in datastru.c
************************************************/
NODE_LIST_T * build_node_list(int type, int n, int **config) {
	NODE_LIST_T *p, *q, *dummy;
	int i;

	p = malloc(sizeof (NODE_LIST_T));
	dummy = p;
	dummy->next = NULL;

	if (type == 0) {
		for (i = 0; i < Node_n; i++) {
			if (config[i][n]) {
				q = malloc(sizeof (NODE_LIST_T));
				q->content = node_p[i];
				p->next = q;
				p = q;
			}
		}
	}

	if (type == 1) {
		for (i = 0; i < Node_n; i++) {
			if (config[n][i]) {
				q = malloc(sizeof (NODE_LIST_T));
				q->content = node_p[i];
				p->next = q;
				p = q;
			}
		}
	}

	p->next = NULL;
	q = dummy;
	dummy = dummy->next;
	free(q);
	return(dummy);
}

/***********************************************
Description:
Sets the parameters in the ANCFIS structure to be modifiable.

fixed = 1 --> fixed parameter;
fixed = 0 --> modifiable parameter

Inputs:
N/A

Outputs:
Returns the number of modifiable parameters.

Called by:
mexFunction in ancfismex.c
************************************************/
int set_parameter_mode() {
	int i, modifiable_p_count = 0;
	PARAMETER_LIST_T *p;
	printf("important check is this \n");
	for (i = 0; i < Node_n; i++) {
		if (node_p[i]->parameter == NULL)
			continue;
		for (p = node_p[i]->parameter; p != NULL; p = p->next) {
			p->fixed = 0;
			modifiable_p_count++;
		}
	}
	printf("Modifiable parameters: %d\n", modifiable_p_count);
	return(modifiable_p_count);
}

