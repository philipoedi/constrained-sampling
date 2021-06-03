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
    const std::size_t d = 3; // dim per sample
    MatrixXd xeig(d,n);
    xeig << 1,1,2,2,3,3;
    double res; // result
    // bandwitdh
    std::vector<double> b{0.1,0.1};
    Vector2d beig(b.data());
    std::cout << xeig*beig << std::endl;
    std::cout << xeig << std::endl;
    std::cout <<beig << std::endl;
    Matrix<double,3,2> de;
    de << 1,1,2,2,3,3;
    std::cout << de << std::endl;
    

    // to predict
    std::vector<double> x_pred{1.5,1.5};
    Vector2d x_predeig(x_pred.data());
    //Map<MatrixXd> x_mat(x.data(),n,d);
    Kernel<d,n> k(beig);
    k.fit(de);
    res = k.evaluate(x_predeig);
    std::cout << res << std::endl;
    assert (round(res*10000)/10000 ==  0.0055);
    
    Kernel<d,n> k2;
    // default constructor test and setBandwidth
    k2.setBandwidth(beig);
    k2.fit(xeig);
    res = k2.evaluate(x_predeig);
    std::cout << res << std::endl;
    assert (round(res*10000)/10000 ==  0.0055);
    
    // find optimal bandwidth
    k2.find_optimal_bandwidth("scott");
    k2.find_optimal_bandwidth("silverman");

    
    std::cout << "Test optimal bandwidth" << std::endl;
    res = k2.evaluate(x_predeig);
    std::cout << res << std::endl;
    
    // Kernel estimator test
    std::vector<Vector2d> x_pred_kde{x_predeig,x_predeig};
    std::vector<double> res_kde;
    res_kde.resize(2);
    KernelEstimator<d,n> kdest;
    kdest.setBandwidth(beig);
    kdest.fit(de);
    kdest.predict(x_pred_kde, res_kde);
    for (const auto& r:res_kde){
        assert (round(r*10000)/10000 == 0.0055);
        std::cout << r << std::endl;
    }
    
    // init with bandwidht est method
    KernelEstimator<d,n> kdest2();
    
    // predit over 2d grid
    std::vector<double> lb{0,0};
    std::vector<double> ub{3,3};
    double step{1};
    kdest.predict(lb, ub, step);
    const std::size_t d2{3};
    KernelEstimator<d2,n> kdest3;
    std::vector<std::vector<double>> da;
    da.resize(d2);
    for (int i=0; i<da.size(); i++)
    {
        da[i].resize(n);
        da[i][0] = 0;
        da[i][1] = 1;
    }

    kdest3.fit(da);

    return 0; 
}
