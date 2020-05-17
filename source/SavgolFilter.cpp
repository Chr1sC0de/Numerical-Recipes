#include <NumericalRecipes/statistics/savgol_coefficients.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>


namespace py = pybind11;
namespace nr = NumericalRecipes;

class SavitzkyGolayFilter {

    private:
        int nl;
        int nr;
        int ld;
        int m;
        nr::types::VectorDouble _coeffs;
 
    SavitzkyGolayFilter(
        const int nl_,
        const int nr_,
        const int ld_=0,
        const int m_=2
    :   
        nl(nl_), nr(nr_), ld(ld_), m(m_), _coeffs(nl+nr+1)
    ){

        int c = nl+nr+1;
        nr::statistics::savgol(_coeffs);
    }

    std::vector<double> __call__(std::vector<double> input){
        std::vector<double> output(input.size() - nl - nr, 0);
        for (int i=0; i<output.size(); i++){
            for (int j=0; j<_coeffs.size(); j++)
                output[i] += _coeffs[j] * input[i+j];
        }
        return output;
    }

    int n_left(){return nl;}
    int n_right(){return nr;}
    int derivative(){return ld;}
    int polynomial(){return m;}
    std::vector<double> coeffs(){
        std::vector<double> output(_coeffs.size());
        for (int i=0; i<_coeffs.size(); i++)
            output[i] = _coeffs[i];
        return output;
    }

};

PYBIND11_MODULE(SavgolFilter, m){

    py::class_<SavitzkyGolayFilter>(m, "SavitzkyGolayFilter")
        .def(py::init<int, int, int, int>(), py::arg("nl"), py::arg("nr"), py::arg("ld"), py::arg("m"))

}

