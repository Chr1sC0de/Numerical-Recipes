#include <NumericalRecipes/base.hpp>

namespace NR = NumericalRecipes;

int main(){

    double carray[3] = {10,10,10};
    double cmatrix[3][3];

    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            cmatrix[i][j] = (double) i*j;

    NR::Vector_t<double> vector_empty;
    NR::Vector_t<double> vector_const_size(10);
    NR::Vector_t<double> vector_const_value(10, 100.0);
    NR::Vector_t<double> vector_const_cval(10, carray[0]);
    NR::Vector_t<double> vector_from_rhs(vector_const_cval);
    NR::Vector_t<double> vector_assigned = vector_from_rhs;

    NR::Matrix2d_t<double> matrix_empty;
    NR::Matrix2d_t<double> matrix_empty_size(3,3);
    NR::Matrix2d_t<double> matrix_constant_value(3,3,10);
    NR::Matrix2d_t<double> matrix_constant_cval(3,3, cmatrix[0]);
    NR::Matrix2d_t<double> matrix_from_rhs(matrix_constant_cval);
    NR::Matrix2d_t<double> matrix_assigned = matrix_from_rhs;

    return 0;

}