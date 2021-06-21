#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <Eigen/Dense>
#include <vector>
#include "utils.hpp"

using namespace Eigen;


enum ConstraintClass {Quadratic, Linear};
enum ConstraintType {Equality, Inequality};


/**
    Auxilliary data used to specifie constraints, such as coefficients,
    constype("ineq" or "eq") and type "linear" or "quadratic"
    @tparam n Dimension N of vector
*/
template<std::size_t n>
struct ConstraintCoeffs {
    /// 
    Matrix<double, n, 1> coeffs = Matrix<double, n, 1>::Zero();
    Matrix<double, n, 1> q = Matrix<double, n, 1>::Zero();
    Matrix<double, n, n> P = Matrix<double, n, n>::Zero();
    double cons;
    double sign;
    double r{0};
    std::string type;
    std::string constype;
};

/**
    Evaluates a linear constraint, function can be used with nlopt
    @param x Current location at which to evaluate constraint 
    @param grad Vector of current gradient
    @param data Auxilliary data to be used in calculations
    @return Constraint value
*/
template<std::size_t n> 
double linearConstraint(const std::vector<double>& x, std::vector<double> &grad, void*data)
{
    typedef Matrix<double, n, 1> vec;
    ConstraintCoeffs<n> *c = (ConstraintCoeffs<n>*) data;
    vec x_vec(x.data());
    if (!grad.empty()){
        utils::copyEig2Vec(c->coeffs, grad);
    }
    std::cout << "slack vals: \n" << x_vec.transpose() * c->coeffs - c->cons << std::endl;
    return x_vec.transpose() * c->coeffs - c->cons;
}

/**
    Evaluates a linear constraint, function can be used with nlopt
    @param x Current location at which to evaluate constraint 
    @param grad Vector of current gradient
    @param data Auxilliary data to be used in calculations
    @return Constraint value
*/
template<std::size_t n> 
double linearConstraint(const Matrix<double,n,1> &x, Matrix<double,n,1> &grad, void*data)
{
    typedef Matrix<double, n, 1> vec;
    ConstraintCoeffs<n> *c = (ConstraintCoeffs<n>*) data;
    if (!grad.empty()){
        grad = c->coeffs;
        //utils::copyEig2Vec(c->coeffs, grad);
    }
    return x.transpose() * c->coeffs - c->cons;
}

/**
    Evaluates a linear constraint, function can be used with nlopt
    @param x Current location at which to evaluate constraint 
    @param grad Vector of current gradient
    @param data Auxilliary data to be used in calculations
    @return Constraint value
*/
template<std::size_t n> 
double linearConstraint(const Matrix<double,n,1> &x, void*data)
{
    ConstraintCoeffs<n> *c = (ConstraintCoeffs<n>*) data;
    return x.transpose() * c->coeffs - c->cons;
}

/**
    Evaluates a quadratic constraint, function can be used with nlopt
    @param x Current location at which to evaluate constraint 
    @param grad Vector of current gradient
    @param data Auxilliary data to be used in calculations, of type ConstraintCoeffs that
        specifies coefficients of constraint
    @return Constraint value
*/
template<std::size_t n> 
double quadraticConstraint(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    typedef Matrix<double, n, 1> vec;
    double res = 0;
    ConstraintCoeffs<n> *c = (ConstraintCoeffs<n>*) data;
    vec x_vec(x.data());
    // 0.5 * x.T@P@x + q.T@x+ r
    if (!grad.empty()){
        utils::copyEig2Vec(c->P.transpose()*x_vec + c->q, grad);
    }
    res += 0.5 * x_vec.transpose() * c->P * x_vec; 
    res += x_vec.transpose()*c->q;  
    std::cout << "x: " << x_vec << std::endl;
    std::cout << "result: " << res -  c->r << "r: "<< c->r << std::endl;
    return res - c->r;
}

template<std::size_t n>
double quadraticConstraint(const Matrix<double,n,1> &x, void*data)
{
    ConstraintCoeffs<n> *c = (ConstraintCoeffs<n>*) data;
    double res{0};
    res += 0.5 * x.transpose() * c->P * x;
    res += x.transpose()*c->q;

    return res - c->r;
}
/*
template<std::size_t n>
class TangentSpace {

    public:
        
        
    private:


};

struct tangentSpaceFinder
{
    Matrix<double,n,m> jac;


};


void tangentSpaceConstraint(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data)
{

    if grad{
        

    }

    Matrix<double,n,m> theta(x[1:].data());

    Matrix<double,m,n> jac = jacobian(x);
    
    Matrix<double,m,m> res1;
    res1 = jac*theta;
    res2 = theta.T@theta - identity<m,m>;
    res1.flatten();
    res2.flatten();
    res = [res1,res2];
    return res;



    //

}
*/

#endif
