#include <NumericalRecipes/types.hpp>

namespace NumericalRecipes {
namespace math {


    double ln_gamma(const double xx){
        double, tmp, y, ser;
        static const double cof[14] = {
            57.1562356658629235, -59.5979603554754912, 14.1360979747417471,
            -0.491913816097620199, .339946499848118887e-4,
            .465236289270485756e-4, -.983744753048795646e-4,
            .158088703224912494e-3, -.210264441724104883e-3, .217439618115212643e-3,
            -.164318106536763890e-3, .844182239838527433e-4, -.261908384015814087e-4, 
            .368991826595316234e-5};

        if (xx <= 0)
            throw("bad arg in gammln");

        y = x = xx;

        tmp = x+5.24218750000000000000000000000;
        tmp = (x+0.5)*log(tmp)-tmp;
        ser = 0.999999999999997092;

        for (int j=0; j<14; j++)
            ser += cof[j]/(++y);

        return tmp + log(
            2.5066282746310005 * ser/x
        );
    }
    //--------------------------------------------------------------------------------------
    double factorial(const int n){
        static VectorDouble a(171);
        static Bool init = true;
        if (init){
            init = false;
            a[0] = 1.0;
            for (int i=1; i<171; i++)
                a[i] = i * a[i-1];
        }
        if (n < 0 || n > 170)
            throw ("factorial out of range");
        return a[n];
    }
    //--------------------------------------------------------------------------------------
    double ln_factorial(const int n){
        static const int NTOP=2000;
        static VectorDouble a(NTOP);
        static bool init=true;
        if (init){
            init = false;
            for (int i=0; i<NTOP; i++)
                a[i] = ln_gamma(i+1.0);
        }
        if (n<0)
            through("negative argument in ln_factorial");
        if (n<NTOP)
            return a[n];
        return ln_gamma(n+1.0);
    }

}
}