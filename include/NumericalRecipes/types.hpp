#pragma once
#include  <NumericalRecipes/base.hpp>

namespace NumericalRecipes {
namespace types {

    using namespace std;

    typedef complex<int> ComplexInt;
    typedef complex<float> ComplexFloat;
    typedef complex<double> ComplexDouble;
    //----------------------------------------------------------------------------------------------------------------------
    typedef const Vector_t<int> VectorInt_I;
    typedef Vector_t<int> VectorInt, VectorInt_O, VectorInt_IO;

    typedef const Vector_t<float> VectorFloat_I;
    typedef Vector_t<float> VectorFloat, VectorFloat_O, VectorFloat_IO;

    typedef const Vector_t<double> VectorDouble_I;
    typedef Vector_t<double> VectorDouble, VectorDouble_O, VectorDouble_IO;

    typedef const Vector_t<ComplexInt> VectorComplexInt_I;
    typedef Vector_t<ComplexInt> VectorComplexInt, VectorComplexInt_O, VectorComplexInt_IO;

    typedef const Vector_t<ComplexFloat> VectorComplexFloat_I;
    typedef Vector_t<ComplexFloat> VectorComplexFloat, VectorComplexFloat_O, VectorComplexFloat_IO;

    typedef const Vector_t<ComplexDouble> VectorComplexDouble_I;
    typedef Vector_t<ComplexDouble> VectorComplexDouble, VectorComplexDouble_O, VectorComplexDouble_IO;
    //----------------------------------------------------------------------------------------------------------------------
    typedef const Matrix2d_t<int> Matrix2dInt_I;
    typedef Matrix2d_t<int> Matrix2dInt, Matrix2dInt_O, Matrix2dInt_IO;

    typedef const Matrix2d_t<float> Matrix2dFloat_I;
    typedef Matrix2d_t<float> Matrix2dFloat, Matrix2dFloat_O, Matrix2dFloat_IO;

    typedef const Matrix2d_t<double> Matrix2dDouble_I;
    typedef Matrix2d_t<double> Matrix2dDouble, Matrix2dDouble_O, Matrix2dDouble_IO;

    typedef const Matrix2d_t<ComplexInt> Matrix2dComplexInt_I;
    typedef Matrix2d_t<ComplexInt> Matrix2dComplexInt, Matrix2dComplexInt_O, Matrix2dComplexInt_IO;

    typedef const Matrix2d_t<ComplexFloat> Matrix2dComplexFloat_I;
    typedef Matrix2d_t<ComplexFloat> Matrix2dComplexFloat, Matrix2dComplexFloat_O, Matrix2dComplexFloat_IO;

    typedef const Matrix2d_t<ComplexDouble> Matrix2dComplexDouble_I;
    typedef Matrix2d_t<ComplexDouble> Matrix2dComplexDouble, Matrix2dComplexDouble_O, Matrix2dComplexDouble_IO;

}
}