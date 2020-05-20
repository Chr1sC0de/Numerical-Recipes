#include <NumericalRecipes/statistics/savgol_coefficients.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>


namespace py = pybind11;
namespace recipes = NumericalRecipes;

class SavitzkyGolayFilter {

    private:

        int nl;
        int nr;
        int ld;
        int m;
        recipes::types::VectorDouble _coeffs;
    
    public:

        SavitzkyGolayFilter(
            const int nl_,
            const int nr_,
            const int ld_=0,
            const int m_=2
        ):
            nl(nl_), nr(nr_), ld(ld_), m(m_), _coeffs(nl+nr+1)
        {
            int np = nl+nr+1;
            recipes::statistics::savgol(_coeffs, np, nl, nr, ld, m);
        }

        void __call__(double * ptr, double * output, int size){

            int output_size(size - nl - nr);
            int coeff_size(_coeffs.size());

            for (int i=0; i<output_size; i++){
                for (int j=0; j<coeff_size; j++)
                    output[i] += _coeffs[j] * ptr[i+j];
            }
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
        .def(
            py::init<int, int, int, int>(),
            py::arg("nl"), py::arg("nr"), 
            py::arg("ld"), py::arg("m")
        ) 
        .def(
            "__call__", 
            [](SavitzkyGolayFilter & self, py::array_t<double> input){
                py::buffer_info buffer = input.request();
                if (buffer.ndim != 1)
                    throw  "input must be a vector";
                
                int size = buffer.size;

                double * ptr = (double *) buffer.ptr;

                int output_size(size - self.n_left() - self.n_right());

                double * output_vector = new double[output_size];

                for (int i=0; i<output_size; i++)
                    output_vector[i] = 0.0;

                self.__call__(ptr, output_vector, size);

                py::capsule free_when_done(output_vector, [](void *f){
                        double *output_vector = reinterpret_cast<double *>(f);
                        delete [] output_vector ;
                });

                return py::array_t<double> (
                    {output_size},
                    {sizeof(double)},
                    output_vector,
                    free_when_done
                );
            }
        )
    ;

}

