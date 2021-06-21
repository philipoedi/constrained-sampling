#ifndef OBJECTIVES_H
#define OBJECTIVES_H
#include <iostream>
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;



/**
    Specifies the vector of Biases in Biased optimization.
    @tparam n Dimension N of vector
*/
template<std::size_t n>
struct Bias {
    /// Eigen Vector representing the bias
    Matrix<double, n, 1> x0;
};



template<std::size_t n, std::size_t m, std::size_t l>
struct SlackData {
    bool state{false};
    Matrix<double, n+m+l+l, 1> a = Matrix<double, n+m+l+l,1>::Zero();
    Matrix<double, n+m+l+l, 1> g = Matrix<double, n+m+l+l,1>::Zero();
};

/**
    Initializes slack variables to use for slack optimization
    @param sl SlackData struct containing vectors of decision variables
*/
template<std::size_t n, std::size_t m, std::size_t l>
void initSlackData(SlackData<n,m,l>& sl)
{
    for (int i=0; i<n+m+l+l; i++)
    {
        if (i >= n && i <n+m)
        {
            sl.a(i) = 1;
            sl.g(i) = 1;
        }
        else if (i >= n+m+l)
        {
            sl.a(i) = 1;
            sl.g(i) = 1;
        };
    };
};

/**
    Evaluates objective function in Biased optimization
    @param x Current location at which to evaluate objective
    @param grad Vector of current gradient
    @param data Auxilliary data to be used in calculations, of type Bias
    @return Objective value
*/
template<std::size_t n> 
double BiasedObjective(const std::vector<double>& x, std::vector<double>& grad, void*data )
{
    typedef Matrix<double, n, 1> vec;
    Bias<n> *b = (Bias<n>*) data;
    vec x_vec(x.data());
    vec x_x0 = x_vec - b->x0;
    if(!grad.empty()){       
        utils::copyEig2Vec(2*x_x0, grad);
/*    for (std::size_t i = 0; i<x_x0.size() ;++i){
        grad[i] = 2*x_x0[i];
    }*/
    }   
    std::cout << "oobj\n" << x_vec<<std::endl;

    std::cout << x_x0.transpose()*x_x0 << std::endl;    
    return x_x0.transpose()*x_x0;    
}

/**
    Evaluates objective function in slack optimization
    @param x Current location at which to evaluate objective
    @param grad Vector of current gradient
    @param data Auxilliary data to be used in calculations
    @return Objective value
*/
template<std::size_t n, std::size_t m, std::size_t l> 
double slackObjective(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    typedef Matrix<double, n+m+l+l,1> vec;
    SlackData<n,m,l>* u = (SlackData<n,m,l>*) data;
    vec x_vec(x.data());
    if(!grad.empty()){
        utils::copyEig2Vec(u->a, grad);
    }
    return u->a.transpose()*x_vec;
}
#endif

