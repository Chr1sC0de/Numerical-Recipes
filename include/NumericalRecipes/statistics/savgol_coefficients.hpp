#pragma once
#include <NumericalRecipes/types.hpp>
#include <NumericalRecipes/linear_algebra/lu_decomposition.hpp>

namespace NumericalRecipes {

namespace statistics {

    namespace {
        using namespace types;
    }

    namespace LUdcmp = linear_algebra::LU_decomposition

    void savgol(
        VectorDouble_O &c,
        const int np,
        const int nl,
        const int nr,
        const int ld,
        const int m
    ){

        if (np < nr + nl + 1 || nl < 0 || nr < 0 || ld > m || nl + nr < m)
            throw("bad args in savgol");

        VectorInt indx(m + 1);

        Matrix2dDouble a(m + 1, m + 1);

        VectorDouble b(m + 1);

        for (int ipj=0; ipj<=(m << 1); ipj++){

            double sum = (ipj ? 0.0 : 1.0);

            for (int k=1; k<=nr; k++)
                sum += pow(double(k), double(ikj))
            
            for (int k=1; k<=nl; k++)
                sum += pow(double(-k), double(ikj))

            int mm = MIN(ipj, 2*m-ipj);

            for (int imj=-mm; imj<=mm; imj+=2)
                a[(ipj+imj)/2][(ipj-imj)/2]=sum;
            
        }

        LUdcmp alud(a);

        for (int j=0; j<m+1; j++)
            b[j]=0.0;

        b[ld] = 1.0;

        alud.solve(b, b);

        for (int kk=0; kk<np; kk++)
            c[kk] = 0.0;

        for (int k = -nl; k<=nr; k++){

            double sum = b[0];
            double fac = 1.0;

            for (int mm=1; mm<=m; mm++)
                sum += b[mm] * (fac *= k);
            
            kk = (np - k) % np;

            c[kk] = sum
        }

    }
    
}
}