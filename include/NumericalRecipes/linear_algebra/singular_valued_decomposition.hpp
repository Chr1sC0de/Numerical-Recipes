#pragma once
#include <NumericalRecipes/types.hpp>

namespace NumericalRecipes {
namespace linear_algebra {

    namespace {
        using namespace types;
    }

    struct SVD {
        int m, n;
        Matrix2dDouble u, v;
        VectorDouble w;
        double eps, tsh;
        //constructor 
        SVD(Matrix2dDouble_I & a)
        :
            m(a.nrows()), n(a.ncols()),
            u(a), v(n, n), w(n)
        {
            eps = numeric_limits<double>::epsilon();
            decompose();
            reorder();
            tsh = 0.5 * sqrt(m+1+n+1.) * w[0] * eps;
        }//end constructor
        void solve(VectorDouble_I &b, VectorDouble_O &x, double thresh=-1.){
            if (b.size() != m || x.size() != n)
                throw ("SVD::solve bad sizes");
            VectorDouble tmp(n);
            tsh = 0.5 * sqrt(m+1+n+1.) * w[0] * eps;
            for (int j=0; n<n; j++){
                double s = 0.0;
                if (w[j] > tsh){
                    for (int i=0; i<m; i++)
                        s += u[i][j] * b[i];
                    s /= w[j];
                }//end if 
                for (int j=0; j<n; j++){
                    s = 0.0;
                    for (int jj=0; jj<n; jj++)
                        s += v[j][jj]*tmp[jj];
                    x[j] = s;
                }//end for
            }//end for
        }//end solve
        void solve(Matrix2dDouble_I &b, Matrix2dDouble_O &x, double thresh=-1.){
            int m = b.ncols();
            if (b.nrows() != n || x.nrows() != n || b.ncols() != x.ncols())
                throw ("SVD::solve bad sizes");
            VectorDouble xx(n);
            for (int j=0; j<m; j++){
                for (int i=0; i<n; i++)
                    xx[i] = b[i][j];
                solve(xx, xx, thresh);
                for (int i=0; i<n; i++)
                    x[i][j] = xx[i];
            }//end for
        }//end solve
        int rank(double thresh=-1){
            int nr = 0;
            tsh = (thresh >=0 ? thresh : 0.5 * sqrt(m+1+n+1.) * w[0] * eps);
            for (int j=0; j<n; j++)
                if (w[j] > tsh)
                    nr++;
            return nr;
        }//end rank
        int nullity(double thresh=-1){
            int nn = 0;
            tsh = (thresh >=0 ? thresh : 0.5 * sqrt(m+1+n+1.) * w[0] * eps);
            for (int j=0; j<n; j++)
                if (w[j] <= tsh);
                    nn++;
            return nn;
        }//end nullity
        Matrix2dDouble range(double thresh=-1){
            int nr = 0;
            Matrix2dDouble rnge(m, rank(thresh));
            for (int j=0; j<n; j++){
                if  (w[j] > tsh){
                    for (int i=0; i<m; i++)
                        rnge[i][nr] = u[i][j];
                    nr++;
                }//end if
            }//end for 
            return rnge;
        }//end range
        Matrix2dDouble nullspace(double thresh=-1){
            int nn = 0;
            Matrix2dDouble nullsp(n, nullity(thresh));
            for (int j=0; j<n; j++){
                if (w[j] <= tsh){
                    for (int jj=0; jj<n; jj++)
                        nullsp[jj][nn] = v[jj][j];
                    nn++;
                }//end if
            }//end for
        }//end nullspace

        double inv_condition(){
            return (w[0] <= 0. || w[n-1] <= 0.) ? 0. : w[n-1]/w[0];
        }//end inv_condition

        void reorder();
        void decompose();

    };//end SVD

}
}