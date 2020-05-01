#include <NumericalRecipes/types.hpp>

namespace NumericalRecipes {
namespace statistics {
    namespace {
        using namespace types;
    }
    void savgol(
        VectorDouble_O &c, const int nl, const int nr,
        const int ld, const int m
    ){
        double fac, sum;
        if (np < nr + nl + 1 || nl < 0 || nr < 0 || ld > m || nl + nr < m)
            throw("bad args in savgol");
        VectorInt indx(m + 1);
        Matrix2dDouble a(m + 1, m + 1);
    }
    
}
}