/* Automatically generated function. */
#include <math.h>

int n_in = 1;
int n_out = 1;
int in_nrow[] = {1};
int in_ncol_[] = {1};
int out_nrow[] = {1};
int out_ncol_[] = {1};
int init(int *n_in_, int *n_out_){
*n_in_ = n_in, *n_out_ = n_out;
return 0;
}

int getInputSize(int n_in, int *n_row, int *n_col){
*n_row = in_nrow[n_in]; *n_col = in_ncol_[n_in];
return 0;
}

int getOutputSize(int n_out, int *n_row, int *n_col){
*n_row = out_nrow[n_out]; *n_col = out_ncol_[n_out];
return 0;
}

int evaluate(const double** x, double** r){
double i_0=x[0][0];
double i_2=(2*i_0);
i_3=(i_0+i_2);
i_5=(3*i_3);
double i_6=(i_3+i_5);
r[0][0]=i_6;
return 0;
}

