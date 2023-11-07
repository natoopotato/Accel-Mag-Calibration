#ifndef ACCEL_CALC_MATHLIB_H
#define ACCEL_CALC_MATHLIB_H

void Multiply_Matrices(double *C, double *A, int nrows, int ncols,
                       double *B, int mcols);
void Matrix_x_Its_Transpose(double *C, double *A, int nrows, int ncols);
void Get_Submatrix(double *S, int mrows, int mcols,
                   double *A, int ncols, int row, int col);
int Choleski_LU_Decomposition(double *A, int n);
int Choleski_LU_Inverse(double *LU,  int n);
void Interchange_Columns(double *A, int col1, int col2, int nrows, int ncols);
void Interchange_Rows(double *A, int row1, int row2, int ncols);

int Hessenberg_Form_Elementary(double *A, double* S, int n);
int QR_Hessenberg_Matrix(double *H, double *Q, double *eigen_real,
                         double *eigen_imag, int n, int iterations);
void Transpose_Square_Matrix( double *A, int n );


#endif //ACCEL_CALC_MATHLIB_H
