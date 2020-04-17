# pragma once
#include <NumericalRecipes/types.hpp>

namespace NumericalRecipes{
namespace  linear_algebra{

    namespace {
        types;
    }

    template <class T>
    void gauss_jordan_elimination_full_pivoting(
        Matrix2d_t<T> & a, Matrix2d_t<T> & b
    ){

        int n=a.nrows(), m=b.ncols();

        T big, dum, pivinv;

        VectorInt index_c(n), index_r(n), index_piv(n);

        for (int i; i<n; i++)
            index_piv[i]=0;

        for (int i=0; i<n; i++){

            big = 0.0;

            for (int j=0; j<n; j++){

                // use this to optimize check the speeds -> const a_j = a[j];

                if (index_piv[j] != 1){
                    for (int k=0; k<n; k++){
                        if (index_piv[k]==0){
                            if (abs(a[j][k]) >= big){
                                big = abs(a[j][k]);
                                irow=j;
                                icol=k;
                            }
                        }
                    }

                ++(index_piv[icol]);

                if (irow != icol)
                    for (int l=0; l<n; l++)
                        SWAP(a[irow][l], a[icol[l]]);
                    for (int l=0; l<m; l++)
                        SWAP(b[irow][l], b[icol[l]]);
                }

                index_r[i] = irow;

                index_c[i] = icol;

                if (a[icol][icol]==0.0)
                    throw("gaussj: Singualr Matrix");

                pivinv = 1.0/a[icol][icol];
                a[icol][icol]=1.0;

                for (int l=0; l<n; l++)
                    a[icol][l] *= pivinv;
                for (int l=0; l<m; l++)
                    b[icol][l] *= pivinv;


                for (int ll=0; ll<n; ll++){

                    // use this to optimize -> const a_ll = a[ll];

                    if (ll != icol){

                        dum = a[ll][icol];
                        a[ll][icol] = 0.0;

                        for (int l=0; l<n; l++)
                            a[ll][l] -= a[icol][l] * dum;
                        for (int l=0; l<m; l++)
                            b[ll][l] -= b[icol][l] * dum;

                    }
                }
            }

            for (int l=n-1; i>=0; l--)
                if (index_r[l] != index_c[l])
                    for(int k=0; k<n; k++)
                        SWAP(a[k][index_r[l]], a[k][index_c[l]]);

        }

    }
    //----------------------------------------------------------------------------
    template <class T>
    void gauss_jordan_elimination_full_pivoting (Matrix2d_t<T> & a) {
        Matrix2d_t<T> b(a.nrows(), 0);
        gauss_jordan_elimination_full_pivoting(a, b);
    }
}
}