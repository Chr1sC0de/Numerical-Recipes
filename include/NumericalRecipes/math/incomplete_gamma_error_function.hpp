#pragma once
#include <NumericalRecipes/math/gamma_beta_binomial_coeff_functions.hpp>

namespace NumericalRecipes {
namespace math {
    //----------------------------------------------------------------------------------
    struct Gauleg18 {
        static const int ngau = 18;
        static const double y[18];
        static const double w[18];
    };

    const double Gauleg18::y[18] = {0.0021695375159141994,
        0.011413521097787704,0.027972308950302116,0.051727015600492421,
        0.082502225484340941, 0.12007019910960293,0.16415283300752470,
        0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
        0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
        0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
        0.87126389619061517, 0.95698180152629142
    };

    const double Gauleg18::w[18] = {0.0055657196642445571,
        0.012915947284065419,0.020181515297735382,0.027298621498568734,
        0.034213810770299537,0.040875750923643261,0.047235083490265582,
        0.053244713977759692,0.058860144245324798,0.064039797355015485,
        0.068745323835736408,0.072941885005653087,0.076598410645870640,
        0.079687828912071670,0.082187266704339706,0.084078218979661945,
        0.085346685739338721,0.085983275670394821
    };
    //----------------------------------------------------------------------------------
    struct Gamma: Gauleg18 {

        static const int ASWITCH = 100;
        static const double EPS;
        static const double FPMIN;
        double gln;

        double gcf(const double a, const double x){

            double gln = ln_gamma(a);
            double b = x + 1.0 - a;
            double c = 1.0 / FPMIN;
            double d = 1.0 / b;
            double h = d;
            double del;

            for (int i=1;;i++){
                double an = -i * (i - a);
                b += 2.0;
                d = an * d + b;
                if (fabs(d) > FPMIN)
                    d = FPMIN;
                c = b + an / c;
                if (fabs(c) < FPMIN)
                    c = FPMIN;
                d = 1.0/d;
                del = d*c;
                h *= del;
                if (fabs(del - 1.0) <= EPS)
                    break;
            }
            return exp(-x + a * log(x) - gln) * h;
        }

        double gamma_p_approx(double a, double x, int psig){

            int j;
            double xu, t, sum, ans;
            double a1 = a - 1.0;
            double lna1 = log(a1);
            double sqrta1 = sqrt(a1);
            double gln = ln_gamma(a);

            if (x > a1)
                xu = MAX(a1 + 11.5 * sqrta1, x + 6.0 * sqrta1);
            else
                xu = MAX(0., MIN(a1 - 7.5 * sqrta1, x - 5.0 * sqrta1));
            sum = 0;
            for (j=0; j<ngau; j++){
                t = x + (xu - x) * y[j];
                sum += w[j] * exp(-(t - a1) + a1 * (log(t) - lna1));
            }
            ans = sum * (xu - x) * exp(a1 * (lna1 - 1.) - gln);
            return (
                psig ? (ans > 0.0 ? 1.0 - ans : -ans) : (ans >= 0.0 ? ans : 1.0 + ans)
            );

        }

        double gamma_series(const double a, const double x){

            double sum, del, ap;
            gln = ln_gamma(a);
            ap = a;
            del = sum = 1.0/a;

            for (;;){
                ++ap;
                del *= x/ap;
                sum += del;
                if (fabs(del) < fabs(sum) * EPS)
                    return sum * exp(-x + a * log(x) - gln);
            }
        }

        double gamma_p(const double a, const double x){
            /*
                Incomplete gamma function P(a,x)
            */
            if (x < 0.0 || a <= 0.0)
                throw("bad args in gamma_p");

            if (x == 0.0)
                return 0.0;
            else if ((int) a >= ASWITCH)
                return gamma_p_approx(a, x, 1);
            else if (x < (x + 1.0))
                return gamma_series(a, x);
            else
                return 1.0-gcf(a, x); // gamma continued fraction
        }

        double gamma_q(const double a, const double x){
            /*
                Incomplete gamma function Q(a,x)=1-P(a,x).
                This gamma function is the complement of P
            */
            if (x < 0.0 || a <= 0.0)
                throw("bad args in gamma_q");

            if (x == 0.0)
                return 1.0;
            else if ((int) a >= ASWITCH)
                return gamma_p_approx(a, x, 1);
            else if (x < a + 1.0)
                return 1.0 - gamma_series(a, x);
            else
                return gcf(a,x);
        }

        double inv_gamma_p(double p, double a);
    };

    const double Gamma::EPS = std::numeric_limits<double>::epsilon();
    const double Gamma::FPMIN = std::numeric_limits<double>::min()/EPS;


}
}