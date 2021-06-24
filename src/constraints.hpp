#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <Eigen/Dense>
#include <vector>
#include "utils.hpp"
#include <algorithm>

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

template<std::size_t n>
double evaluateConstraint(std::vector<double> &x, ConstraintCoeffs<n> &con){
    double slack{0};
    std::vector<double> grad;
    if (con.constype == "linear"){
        slack = linearConstraint<n>(x, grad, &con);
    } else if (con.constype == "quadratic"){
        slack = quadraticConstraint<n>(x, grad, &con);
    } else {
        std::cout << "wrong constraint type used. set .cosntype to either 'quadratic' or 'linear'" << std::endl;
        exit(1);
    }
    return slack;
}


template<std::size_t n>
bool isFeasible(std::vector<double> &x, ConstraintCoeffs<n> &cons){
   double slack = evaluateConstraint<n>(x, cons);
   if (cons.type == "eq") {
        return (slack ==  0) ? true : false;
    } else if (cons.type == "ineq"){
        return (slack <= 0) ? true : false; 
    } else {
        std::cout << "Choose ConstraintCoeffs.type = {'eq','ineq'} " << std::endl;
        exit(1);
    }
}


template<std::size_t n>
bool isFeasiblePtr(std::vector<double> &x, ConstraintCoeffs<n> *cons){
   double slack = evaluateConstraint<n>(x, *cons);
   if (cons->type == "eq") {
        return (slack ==  0) ? true : false;
    } else if (cons->type == "ineq"){
        return (slack <= 0) ? true : false; 
    } else {
        std::cout << "Choose ConstraintCoeffs.type = {'eq','ineq'} " << std::endl;
        exit(1);
    }
}


template<std::size_t n>
bool isFeasibleM(std::vector<double> &x, std::vector<ConstraintCoeffs<n>> &cons){
    if (std::all_of(cons.begin(), cons.end(), std::bind(isFeasible<n>,x, std::placeholders::_1))){
        return true;
    } else {
        return false;
    }
}


template<std::size_t n>
bool isFeasibleM(std::vector<double> &x, std::vector<ConstraintCoeffs<n>*> &cons){
    if (std::all_of(cons.begin(), cons.end(), std::bind(isFeasiblePtr<n>,x, std::placeholders::_1))){
        return true;
    } else {
        return false;
    }
}

template<std::size_t n>
bool boundsCheck(Matrix<double,n,1> &x, Matrix<double,n,1> &lb, Matrix<double,n,1> &ub){
    if ((ub.array()>=x.array()).all() && ((lb.array()<=x.array()).all())){
        return true;
    } else {
        return false;
    }
}

template<std::size_t n>
bool boundsCheckVec(std::vector<double> &x, Matrix<double,n,1> &lb, Matrix<double,n,1> &ub ){
    Matrix<double,n,1> x_eig(x.data());
    return boundsCheck<n>(x_eig,lb,ub);
}


std::vector<double> numericGradient(const std::vector<double> &x, std::function<double(std::vector<double>&)> func, double h){
     std::vector<double> x_temp = x;
     std::vector<double> grad(x.size());
     double f_big, f_small;
     for (int i=0; i<x.size(); i++){
        x_temp[i] += h;
        f_big = func(x_temp);
        x_temp[i] -= 2*h;
        f_small = func(x_temp);
        x_temp[i] += h;
        grad[i] = (f_big - f_small) / (2*h);
     }
     return grad;
}


template<std::size_t n, std::size_t m>
Matrix<double,m,n> numericJacobian(const std::vector<double> &x, std::vector<std::function<double(std::vector<double>&)>> funcs, double h){
    assert (m>0);
    if (m==0) exit(1);
    Matrix<double,m,n> jac;
    for (int i=0; i<funcs.size(); i++){
       jac.row(i) = Matrix<double,1,n>(numericGradient(x,funcs[i],h).data()); 
    }
    return jac;
}

template<std::size_t n>
ConstraintCoeffs<n> createSphere(double radius){
    ConstraintCoeffs<n> sphere;
    sphere.constype = "quadratic"; 
    sphere.type = "eq";
    sphere.r = radius;
    sphere.P = Matrix<double,n,n>::Identity();
    sphere.P *= 2;
    return sphere; 
};


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
