
//#include "data_formatter.c"
//#include "ancfismex.h"  I have done this.. Omkar has done this
#include "ancfis.h"
NODE_T **node_p;
int In_n, Out_n, In_vect_n, Mf_n, Rule_n, Node_n;

///////////////////////////////////////////////////////////////////////////////////////////////////

//void main()
int main()
{
	int n_ip = 1,n_op = 1,n_trn_data = 1,n_chk_data = 0,ip_len =16,ss = 0.001,stepsize=0,n_mf =3,n_epochs = 100,ss_d = 0.8,ss_i = 1.1,debug = 0,random_iterations = 1;

int hhh[]={1, 1, 1, 1, 1};
int hhh2[]={1, 1, 1, 1, 1};
int min_trn_err = 9999;
int min_chk_err = 9999;
int min_trn_err_epoch = 1;
int min_chk_err_epoch = 1;
int avg = 0;
int min_avg = 9999;
int st_trn[] = {}; 
int st_chk[] = {};
int st_trn_epoch[] = {};
int st_chk_epoch[] = {};
int index = 1;

	n_trn_data = data_formatter('nor_synthetic_trn.txt','data.trn',1000,ip_len,1,1,1,1,1,1,1); 
    n_chk_data = data_formatter('nor_synthetic_chk.txt','data.chk',565,ip_len,1,1,1,1,1,1,1);

	a=ancfismex(n_ip,n_op,ip_len,n_mf,n_trn_data,n_chk_data,n_epochs,ss,ss_d,ss_i,debug);
                        
                       
                        st_trn[index] = a(1);
                        st_chk[index] = a(2);
                        st_trn_epoch[index] = a(3);
                        st_chk_epoch[index] = a(4);
                        
                        index = index + 1;

                       if (a(1)<min_trn_err)
			{
                            min_trn_err = a(1);
                            min_trn_err_epoch = a(3)+1;                       
                        //end
			return(1) ;
			}
                        
                        if (a(2)<min_chk_err)
			{	
                        //hhh[] = {ip_len, ss, ss_d, ss_i, n_mf};    
			    	hhh[0] = ip_len;
				hhh[1] = ss;
				hhh[2] = ss_d;
				hhh[3] = ss_i;
				hhh[4] = n_mf;
                            	min_chk_err = a(2);
                            	min_chk_err_epoch = a(4)+1;                       
                        //end
                        
                        avg(random_iterations) = a(2);                        
                       // %                    end
                    return(2) ;
			}

                    
                    if (min_avg > mean(avg))
			{
                        min_avg = mean(avg);
                        //hhh2[] = {ip_len, ss, ss_d, ss_i, n_mf}; 
			hhh2[0] = ip_len;
			hhh2[1] = ss;
			hhh2[2] = ss_d;
			hhh2[3] = ss_i;
			hhh2[4] = n_mf;                        
                    return(0) ;
			}

}
