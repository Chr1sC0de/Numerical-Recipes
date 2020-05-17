#pragma once
#include <NumericalRecipes/types.hpp>

namespace NumericalRecipes {

namespace linear_algebra {

    namespace {
        using namespace types;
    }

    struct LU_decomposition {
        int n;
        // sotres the decomposition
        Matrix2dDouble lu;
        //stores the permutation
        VectorInt index;
        //used by the determinant
        double d;
        //only used by mprove
        Matrix2dDouble_I &aref;
        //constructor
        LU_decomposition(Matrix2dDouble_I & a)
        :
            n(a.nrows()), lu(a), aref(a), index(n)
        {
            const double TINY = 1.0e-40;    
            VectorDouble vv(n);
            d = 1.0;
            int imax;

            for (int i=0; i<n; i++){
                double big=0.0;
                for (int j=0; j<n; j++){
                    double temp = lu[i][j];
                    if (temp > big)
                        big = temp;
                }//end for
                if (big == 0.0)
                    throw("Singular matrix in LUdcmp");
                vv[i] = 1.0/big;
            }//end for

            for (int k=0; k<n; k++){
                double big=0.0;
                for (int i=k; i<n; i++){
                    double temp=vv[i]*abs(lu[i][k]);
                    if (temp > big){
                        big = temp;
                        imax = i;
                    }//end if
                }//end for
                if (k!= imax) {
                    for (int j=0; j<n; j++){
                        double temp = lu[imax][j];
                        lu[imax][j] = lu[k][j];
                        lu[k][j] = temp;
                    }//end for
                    d = -d;
                    vv[imax] = vv[k];
                }//end if
                index[k] = imax;
                if (lu[k][k] == 0.0)
                    lu[k][k] = TINY;
                for (int i=k+1; i<n; i++){
                    double temp = lu[i][k]/=lu[k][k];
                    for (int j=k+1; j<n; j++)
                        lu[i][j] -= temp*lu[k][j];
                }//end for
            }//end for 
        }//end constructor
        void solve(VectorDouble_I &b, VectorDouble_O & x){
            /*
                solve A x = b -> L (U x) = b where L and U
                have been decomposed by the construtor.
                x is the solution and will be stored
            */
            int ii=0;
            double sum; 

            if (b.size() != n || x.size() != n)
                throw ("LU_decomposition::solve bad sizes");

            for (int i=0; i<n; i++)
                x[i] = b[i];

            for (int i=0; i<n; i++){
                int ip = index[i];
                sum = x[ip];
                x[ip] = x[i];
                if (ii != 0)
                    for (int j=ii-1; j<i; j++)
                        sum -= lu[i][j]*x[j];
                else if (sum != 0.0)
                    ii += i+1;
                x[i] = sum;
            }//end for 
            for (int i=n-1; i>=0; i--){
                sum = x[i];
                for (int j=i+1; j<n; j++)
                    sum -= lu[i][j] * x[j];
                x[i] = sum/lu[i][i];
            }//end for
        }//end solve
        void solve (Matrix2dDouble_I & b, Matrix2dDouble_O & x){
            int m = b.ncols();
            if (b.nrows() != n || x.nrows() != n || b.ncols() != x.ncols())
                throw("LU_decomposition::solve bad sizes");
            VectorDouble xx(n);
            for (int j=0; j<m; j++){
                for (int i=0; i<n; i++)
                    xx[i] = b[i][j];
                solve(xx, xx);
                for (int i=0; i<n; i++) 
                    x[i][j] = xx[i];
            }//end for 
        }//end solve
        void inverse(Matrix2dDouble_O &ainv){
            ainv.resize(n, n);
            for (int i=0; i<n; i++){
                for (int j=0; j<n; j++)
                    ainv[i][j] = 0;
                ainv[i][i] = 1.0;
            }//end for
            solve(ainv, ainv);
        }//end inverse
        double det(){
            double dd = d;
            for (int i=0; i<n; i++)
                dd *= lu[i][i];
            return dd;
        }//end det
        void mprove(VectorDouble_I &b, VectorDouble_IO &x){
            VectorDouble r(n);
            for (int i=0; i<n; i++){
                long double sdp = -b[i];
                for (int j=0; j<n; j++){
                    sdp += (long double)aref[i][j] * (long double)x[j];
                    r[i] = sdp;
                }//end for
                solve(r, r);
                for (int i=0; i<n; i++)
                    x[i] -= r[i];
            }//end for
        }//end mprove
    };//end LU_decomposition

}
}