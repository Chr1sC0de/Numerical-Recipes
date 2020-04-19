#include <NumericalRecipes/math/gamma_beta_binomial_coeff_functions.hpp>
#include <NumericalRecipes/math/incomplete_gamma_error_function.hpp>

namespace NR = NumericalRecipes;

int main(){

    double result_1 = NR::math::ln_gamma(15);
    double result_2 = NR::math::factorial(10);
    double result_3 = NR::math::ln_factorial(10);
    double result_4 = NR::math::binomial_coefficient(10, 4);

    return 0;
}