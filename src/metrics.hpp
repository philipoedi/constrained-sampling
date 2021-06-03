#ifndef METRICS_H
#define METRICS_H

#include <iostream>
#include <vector>
#include <algorithm>
#include "optimizer.hpp"

double objective_coverage(const std::vector<double> &x, std::vector<double> &grad, void *data){
    if (!grad.empty()){
        std::fill(grad.begin(), grad.end(), 0);    
        grad[0] = -1;    
    }
    return -x[0];
}

/*
template<std::size_t n> double coverage_constraint(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    ConstraintCoeffs<n>* c = (ConstraintCoeffs<n>*) data;
    typedef Matrix<double, n,1> vec; 
    typedef Matrix<double, n+1,1> vec_1; 
    double t = x[0];
    double s = c->sign;
    vec xstar = c->coeffs;
    vec_1 x_vec(x.data());
    vec x_xstar = x_vec.tail(n) - xstar;
    if (!grad.empty()){
       for (std::size_t i = 1; i < x_vec.size(); ++i){
           grad[i] = s *  2 * x_xstar[i-1];
       }
       grad[0] = 1;
    }
    double out = x_xstar.transpose()*x_xstar;
    std::cout << "------------" <<std::endl;
    std::cout << x_vec << std::endl;
    std::cout << "out" << std::endl;
    std::cout << out << std::endl;
    std::cout << "s*out+t" << std::endl;
    std::cout << s*out+t <<std::endl;
    std::cout << "x_star" << std::endl;
    std::cout << xstar << std::endl;
    std::cout << c->coeffs << std::endl;
    std::cout << s <<std::endl;
    return s * out + t;
};

*/

// Kernel llo function to add, where predict other arg index given -> this index is not included in eval
double loo_entropy()
    
    KernelEstimator<n-1,d> ( ) 
    vec llo_new ; // vector left out
    vec llo_old ;
    for ( i in range) 
        llo_old = llo_new;
        llo_new  = Kernel.data[i] 
        
    
#endif



