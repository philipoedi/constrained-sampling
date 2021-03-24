#include "kde.hpp"
#include <vector>
#include <iostream>
#include <cassert>
#include <limits>
#include <Eigen/Dense>
#include <string>

int main(){
    // training data declaration
    const std::size_t n = 2; // num samples
    const std::size_t d = 2; // dim per sample
    MatrixXd xeig(d,n);
    xeig << 1,1,2,2;
    double res; // result
    // bandwitdh
    std::vector<double> b{0.1,0.1};
    Vector2d beig(b.data());
    // to predict
    std::vector<double> x_pred{1.5,1.5};
    Vector2d x_predeig(x_pred.data());
    //Map<MatrixXd> x_mat(x.data(),n,d);
    kernel<d,n> k(beig);
    k.fit(xeig);
    res = k.evaluate(x_predeig);
    std::cout << res << std::endl;
    assert (round(res*10000)/10000 ==  0.0056);
    kernel<d,n> k2;
    // default constructor test and set_bandwidth
    k2.set_bandwidth(beig);
    k2.fit(xeig);
    res = k2.evaluate(x_predeig);
    std::cout << res << std::endl;
    assert (round(res*10000)/10000 ==  0.0056);
    // find optimal bandwidth
    k2.find_optimal_bandwidth("scott");
    res = k2.evaluate(x_predeig);
    std::cout << res << std::endl;
    // kernel estimator test
    std::vector<Vector2d> x_pred_kde{x_predeig,x_predeig};
    std::vector<double> res_kde;
    res_kde.resize(2);
    kernel_estimator<d,n> kdest;
    kdest.set_bandwidth(beig);
    kdest.fit(xeig);
    kdest.predict(x_pred_kde, res_kde);
    for (const auto& r:res_kde){
        assert (round(r*10000)/10000 == 0.0056);
        std::cout << r << std::endl;
    }
}
