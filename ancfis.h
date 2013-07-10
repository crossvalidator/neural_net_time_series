/*
This is the header file for ANCFIS.

There are 4 data structures used in ANCFIS,
1. NODE_T represents one ANCFIS node.
2. NODE_LIST_T are used to define the connections between nodes
3. PARAMETER_LIST_T are used to store parameters in layers 1 and 4
4. COMPLEX_T is a structure for complex numbers
*/

#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//-----------DATA STRUCTURES
typedef struct node_s {
	int index; //global node index within network
	int layer;
	int local_index; //index of the node within a layer
	int parameter_n; //number of parameters
	struct complex_s *value;
	struct complex_s *de_do;
	double *input_vector; //input vector, used for the input nodes, null otherwise
	struct complex_list_s *dE_dSMF;
	int function_index;	//node function index, used in forward pass (now, this is always the same as layer)
	struct node_list_s *fan_in; //list of fan_in nodes
	struct node_list_s *fan_out; //list of fan_out nodes
	struct parameter_list_s *parameter;	//list of parameters
} NODE_T;

typedef struct node_list_s {
	NODE_T *content;
	struct node_list_s *next;
} NODE_LIST_T;

typedef struct parameter_list_s {
	int fixed;
	double content; //the value of the parameter
	struct complex_s *de_dp;
	struct parameter_list_s *next;
} PARAMETER_LIST_T;

typedef struct complex_s {
	double real, imag;
} COMPLEX_T;

typedef struct complex_list_s {
	struct complex_s *content;
	struct complex_list_s *next;
} COMPLEX_LIST_T;
/*
typedef struct OA_rr {
	double minTrnRMSE;
	double minChkRMSE;
	double minTrnRMSEEpoch;
	double minChkRMSEEpoch;
} OA_RR;
*/
//-----------DATA STRUCTURES

//>>>>>>>>>>>GLOBAL VARIABLES
extern In_n;  /* number of input variables, should be 1 */
extern Out_n; //number of outputs
extern In_vect_n; //number of data points in the input vector
extern Mf_n;   /* number of membership functions along each input */
extern Node_n; /* number of total nodes */
extern Rule_n;    /* number of nodes in the 4-th layer */
extern NODE_T **node_p; //an array of nodes, which holds the whole ANCFIS structure
//>>>>>>>>>>>GLOBAL VARIABLES

//===========FUNCTION PROTOTYPES
//memory allocation
char *create_array();
char **create_matrix();
void free_array();
void free_matrix();

//print to stdio and write to files
void print_array();
void print_matrix();
void write_array();
void write_matrix();

//complex functions in complex.c

COMPLEX_T complex();
COMPLEX_T c_add();
COMPLEX_T c_sub();
COMPLEX_T c_mul();
COMPLEX_T c_conjg();
COMPLEX_T c_div();
double c_abs();
COMPLEX_T c_sqrt();
COMPLEX_T c_sqrt();
COMPLEX_T c_mul_scalar();
COMPLEX_T c_add_scalar();
COMPLEX_T c_dot_product();

COMPLEX_LIST_T * build_complex_list(int n);
COMPLEX_T * first_initialized_weights();
//===========FUNCTION PROTOTYPES

//+++++++++++DEFINITIONS
#define TRAIN_ERR_FILE   "error.trn"
#define CHECK_ERR_FILE   "error.chk"
#define INIT_PARA_FILE   "para.ini"
#define FINA_PARA_FILE   "para.fin"
#define TRAIN_DATA_FILE  "data.trn"
#define CHECK_DATA_FILE  "data.chk"
#define STEP_SIZE_FILE   "stepsize"
#define INPUT_FILE   "inputpara.ini"


/* error type definition, see trn_err.c and  */
/* 0 --> RMSE (root mean square error)
   1 --> NDEI (non-dimensional error index)
   2 --> ARV  (average relative variance)
   3 --> APE  (average percentage error)
   4 --> MAE  (mean absolute error)
   */
#define ERROR_TYPE 0 //used in chk_err.c and trn_err.c
#define PI 3.1415926535897932385E0
#define RANDOM_INIT 1 //random parameter initialization = 1, see initpara.c
//#define R_n 8
//+++++++++++DEFINITIONS
