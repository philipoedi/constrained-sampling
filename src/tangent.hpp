#ifndef TANGENT_H
#define TANGENT_H
#include <iostream>
#include "constraints.hpp"
#include "optimizer.hpp"
#include <Eigen/Dense>
#include <nlopt.hpp>
#include "objectives.hpp"

template<std::size_t n, std::size_t m>
class BaseSampler;

template<std::size_t n, std::size_t m>
class UniformSampler;

#include "sampler.hpp"

using namespace Eigen;
using namespace nlopt;
using namespace std::placeholders;

template<std::size_t n, std::size_t m>
struct TangentConstraintData{
    Matrix<double,m*(n-m)+(n-m)*(n-m) , n*(n-m),RowMajor> A = Matrix<double,m*(n-m)+(n-m)*(n-m),n*(n-m)>::Zero();
    Matrix<double,n*(n-m),1> b = Matrix<double,n*(n-m),1>::Zero();
};

template<std::size_t n, std::size_t m> // n dimension ambient space, m equality constraints
class TangentSpace {

    typedef Matrix<double,n,1> AmbientVector;
    typedef Matrix<double,n-m,1> TangVector; // Vector with dimension of the tangent space 
    typedef Matrix<double,m*(n-m)+(n-m)*(n-m),n*(n-m),RowMajor> ConsMatrix;
    typedef Matrix<double,n*(n-m),1> ThetaFlat;
    typedef Matrix<double,n,n-m> ThetaMatrix;

    public:

        TangentSpace();
        TangentSpace(const AmbientVector &x0, const Matrix<double,m,n> &jac);
        void findTangentSpace(const AmbientVector &x0, const Matrix<double,m,n> &jac); // jacobian of implicit manifold def
        void createConstraintData(const Matrix<double,m,n> &jac);
        void updateConstraintData(const ThetaMatrix &theta);
        TangentConstraintData<n,m> getConstraintData();
        AmbientVector toAmbient(const TangVector &u);
        ConsMatrix getGradient(const ConsMatrix &A);
        ConsMatrix getGradient();
        ThetaFlat getSlack(const ConsMatrix &A, const ThetaFlat &theta, const ThetaFlat &b);
        ThetaFlat getSlack(const ThetaFlat &theta);
        ThetaMatrix getTheta();
        TangVector sampleOnTangent();
        AmbientVector sampleInAmbient();
        template<typename T>
        void setSamplerBounds(const T &lb, const T &ub);

    private:

        TangentConstraintData<n,m> data_;
        ThetaMatrix theta_; // map from lower dimensional tangent space to ambient space
        AmbientVector x0_;
        std::size_t k_{n-m};
        //opt opt_{"AUGLAG_EQ",n*(n-m)};
        opt opt_{"LN_COBYLA",n*(n-m)};
        opt local_opt_{"LD_SLSQP",n*(n-m)};
        UniformSampler<n-m,0> uni_;
        //BiasedOptimizer<n*(n-m)> bopt_;
};

template<std::size_t n, std::size_t m>
void tangentSpaceConstraints(unsigned l, double *result, unsigned k, const double *x, double *grad, void* f_data);


template<std::size_t n, std::size_t m>
TangentSpace<n,m>::TangentSpace(){
    double tol{1e-8};
    opt_.set_xtol_rel(tol);
    local_opt_.set_xtol_rel(tol);
    //opt_.set_local_optimizer(local_opt_);
    
}


template<std::size_t n, std::size_t m>
Matrix<double,n,n-m> TangentSpace<n,m>::getTheta(){
    return theta_; 
}   

template<std::size_t n, std::size_t m>
void TangentSpace<n,m>::createConstraintData(const Matrix<double,m,n> &jac){
    for (int i=0; i < (n-m) ; i++){

        data_.A.block(i*m,i*n,m,n) = jac ;
    }
    for (int i=0; i<(n-m); i++ ){   
        data_.b(i*(n-m)+(m*(n-m))+i) = 1;
    };
}

template<std::size_t n, std::size_t m>
void TangentSpace<n,m>::updateConstraintData(const ThetaMatrix &theta){
    for (int i=0; i<(n-m); i++){
        data_.A.block((i*(n-m))+(m*(n-m)),i*n,(n-m),n) = theta.transpose();
    }
}

template<std::size_t n, std::size_t m>
TangentConstraintData<n,m> TangentSpace<n,m>::getConstraintData(){
    return data_;
}

template<std::size_t n, std::size_t m>
Matrix<double,n*(n-m),1> TangentSpace<n,m>::getSlack(const ConsMatrix &A, const ThetaFlat &theta, const ThetaFlat &b){
    return A*theta - b;
}


template<std::size_t n, std::size_t m>
Matrix<double,m*(n-m)+(n-m)*(n-m),n*(n-m),RowMajor> TangentSpace<n,m>::getGradient(const ConsMatrix &A){
    ConsMatrix grad = A;
    for (int i=0; i<(n-m); i++){
        grad.row(i*(n-m) + m*(n-m)+i) *= 2;
    }
    return grad;
}


template<std::size_t n, std::size_t m>
Matrix<double,n*(n-m),1> TangentSpace<n,m>::getSlack(const ThetaFlat &theta){
    return  getSlack(data_.A, theta, data_.b);
}


template<std::size_t n, std::size_t m>
Matrix<double,m*(n-m)+(n-m)*(n-m),n*(n-m),RowMajor> TangentSpace<n,m>::getGradient(){
    return getGradient(data_.A);
}

template<std::size_t n, std::size_t m>
template<typename T>
void TangentSpace<n,m>::setSamplerBounds(const T &lb, const T &ub){
    uni_.setBounds(lb, ub);
}

template<std::size_t n, std::size_t m>
Matrix<double,n-m,1> TangentSpace<n,m>::sampleOnTangent(){
    return uni_.sample();
}


template<std::size_t n, std::size_t m>
Matrix<double,n,1> TangentSpace<n,m>::sampleInAmbient(){
    return toAmbient(uni_.sample());
}

template<std::size_t n, std::size_t m>
void TangentSpace<n,m>::findTangentSpace(const AmbientVector &x0, const Matrix<double,m,n> &jac){
    createConstraintData(jac);
    x0_ = x0;
    std::vector<double> tols(m*(n-m)+(n-m)*(n-m),1e-8);
    
    opt_.add_equality_mconstraint(tangentSpaceConstraints<n,m>, this, tols);
   // opt_.set_maxeval(50);
    // run 
    std::vector<double> x(n*(n-m),-1);
    std::vector<double> lb(n*(n-m),-10);
    std::vector<double> ub(n*(n-m),10);
/*    UniformSampler<n*(n-m)> uni(lb,ub);
    x = utils::copyEig2Vec(uni.sample());
  */  Bias<n*(n-m)> b;
    ThetaFlat b0(x.data());
    b.x0 = b0;
    double minf;
    opt_.set_min_objective(BiasedObjective<n*(n-m)>, &b);  
    result  res = opt_.optimize(x, minf);
    ThetaMatrix theta(x.data());
    theta_ = theta;
}

template<std::size_t n, std::size_t m>
Matrix<double,n,1> TangentSpace<n,m>::toAmbient(const TangVector &u){
    return x0_ + theta_*u;
}


template<std::size_t n, std::size_t m>
void tangentSpaceConstraints(unsigned l, double *result, unsigned k, const double *x, double *grad, void* f_data){
    TangentSpace<n,m>  *tang = reinterpret_cast<TangentSpace<n,m>*>(f_data);
    Matrix<double,n,n-m> theta_m(x);
    Map<const Matrix<double,n*(n-m),1>> theta(theta_m.data(),theta_m.size());
    tang->updateConstraintData(theta_m);
    /**if (grad){
        Matrix<double,m*(n-m)+(n-m)*(n-m),n*(n-m),RowMajor> grad_eig = tang->getGradient();
        Map<const VectorXd> grad_vec(grad_eig.data(),grad_eig.size());
        std::cout <<"gradient \n"<< grad_vec << std::endl;
        utils::copyEig2Arr(grad_vec, grad); 
    }*/
    Matrix<double,n*(n-m),1> result_eig  = tang->getSlack(theta);
    utils::copyEig2Arr(result_eig, result);
}


template<std::size_t n, std::size_t m>
TangentSpace<n,m> tangentSpaceFromConstraints(const std::vector<ConstraintCoeffs<n>*> &cons_ptrs, std::vector<double> x, double h){
    ConstraintCoeffs<n> * c_ptr;
    std::vector<std::function<double(std::vector<double>&)>> funcs;
    for (int i=0; i<cons_ptrs.size(); i++){
        c_ptr = cons_ptrs[i];
        if (c_ptr->type == "eq"){
            auto f = std::bind(evaluateConstraint<n>,_1, *c_ptr);
            funcs.push_back(f);
        }
    }
    Matrix<double,m,n> jac = numericJacobian<n,m>(x,funcs, h);
    TangentSpace<n,m> tang;
    Map<Matrix<double,n,1>> x_eig(x.data(),n); 
    tang.findTangentSpace(x_eig,jac);
    return tang;
}


#endif
