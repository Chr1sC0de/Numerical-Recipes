#include <NumericalRecipes/io/csv.hpp>
#include <NumericalRecipes/plot/psplot.hpp>
#include <limits>


using namespace NumericalRecipes;


typedef std::unordered_map<std::string, types::VectorDouble> DOHLCV;

int main(){

    std::string dohlcv_file = "resources/BinancBTC_4hr.csv";
    std::string save_file = "D:\\Github\\Numerical-Recipes\\resources\\candle_line_plot.ps";

    io::TableReader csv_reader(dohlcv_file);

    DOHLCV data = csv_reader.get_data();

    types::VectorDouble open = data["open"];

    double max(-10000);

    for (int i=0; i<open.size(); i++){
        if (open[i] > max){
            max = open[i];
        }
    }

    types::VectorDouble datetime = data["datetime"];

    int n=1000;

    // now plot the open data
    types::VectorDouble test_x(n), test_y(n);

    for (int i=0; i<n; i++){
        test_x[i] = 5.*i/(n-1.);
        test_y[i] = exp(-0.5 * test_x[i]);
    }

    double max_datetime(datetime[datetime.size() - 1]);

    for (int i=0; i<open.size(); i++){
        double d_temp = datetime[i];
        double o_temp = open[i];

        datetime[i] = d_temp/max_datetime * 100;

        if (o_temp<0){
            o_temp = 0;
            open[i] = o_temp;
        }

    }

    psplot::Page pg(save_file);
    psplot::Axis ax(pg, 100., 500., 100., 500.);

    ax.setlimits(0., 100., 0., max);

    ax.frame();
    ax.autoscales();
    ax.lineplot(datetime, open);

    return 0;
}