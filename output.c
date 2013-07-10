//#include "ancfis.h"
/*
This file contains functions that 'store' data for further viewing or use.
Here, store means 'write to an array' or 'write to a file'.

There are 3 functions:
1. write_parameter
2. record_parameter
3. restore_parameter
*/

/***********************************************
Description:
Write all node parameters to file 'file_name'.

Inputs:
file_name - the file name

Outputs:
A file.

Called from:
mexFunction in ancfismex.c
************************************************/
void write_parameter(char *file_name) {
	FILE *fp;
	PARAMETER_LIST_T *p;
	int i;

	fp = (FILE *)open_file(file_name, "w");

	for(i = 0; i < Node_n; i++) {
		if(node_p[i]->parameter == NULL)
			continue;
		for(p = node_p[i]->parameter; p != NULL; p = p->next)
			fprintf(fp, "%4.10lf ", p->content);
		fprintf(fp, "\n");
	}

	fclose(fp);
}

/***********************************************
Description:
Record the current parameters in an array.  This is used to remember the
parameters that give the best training error.

Inputs:
*parameter_array - must be preallocated.
	This will store the ANCFIS parameters

Outputs:
Copies the parameter values into *parameter_array.

Called from:
mexFunction in ancfismex.c
************************************************/
void record_parameter(double *parameter_array) {
	int i, j = 0;
	PARAMETER_LIST_T *p;

	for(i = 0; i < Node_n; i++) {
		if(node_p[i]->parameter == NULL)
			continue;
		for(p = node_p[i]->parameter; p != NULL; p = p->next)
			parameter_array[j++] = p->content;
	}
}

/***********************************************
Description:
Restore the parameters in parameter_array (created using function
record_parameter) back to ANFIS.

Inputs:
*parameter_array - must be preallocated.
	This holds the stored ANCFIS parameters from record_parameter.

Outputs:
Operates on node_p[i]->parameter

Called from:
mexFunction in ancfismex.c
************************************************/
void restore_parameter(double *parameter_array) {
	int i, j = 0;
	PARAMETER_LIST_T *p;

	for (i = 0; i < Node_n; i++) {
		if (node_p[i]->parameter == NULL)
			continue;
		for (p = node_p[i]->parameter; p != NULL; p = p->next)
			p->content = parameter_array[j++];
	}
}



