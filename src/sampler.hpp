#ifndef SAMPLER_H
#define SAMPLER_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <random>
#include <functional>
#include <cassert>
#include <memory>


template<std::size_t n>
class BaseOptimizer;

template<std::size_t n>
class BiasedOptimizer;

#include "utils.hpp"
#include "constraints.hpp"
#include "optimizer.hpp"

using namespace Eigen;


template<std::size_t n>
class BaseSampler
{
    public:
        typedef Matrix<double, n, 1> Vector;
        enum {NeedsToAlign = (sizeof(Vector)%16)== 0};
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign); 
        BaseSampler(){};
        template <typename T> BaseSampler(const T &lb, const T &ub);
        void setBounds(const Vector& lb, const Vector& ub);
        void setBounds(const std::vector<double> &lb, const std::vector<double> &ub);
        std::pair<Vector,Vector> getBounds();
        virtual void run(int n_iter){};
        virtual void run(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub){};
        virtual void runOnTangent(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub){};
        std::vector<std::vector<double>> results();
        std::vector<std::vector<double>> samples();
        bool checkFeasible(const std::vector<double> &x);
        bool checkFeasible(const Vector &x);
        void addConstraints(std::vector<ConstraintCoeffs<n>> &cons);
        //void setOptimizer(BaseOptimizer<n> * opt);
        void setOptimizer(std::string method, std::vector<double> &lb, std::vector<double> &ub);
        void saveResults(std::string name);
        void saveSamples(std::string name);
        bool hasOptimizer();
        std::vector<double> optimize(std::vector<double> &seed);
        void reset();

    protected:
        Vector lb_, ub_, x_;
        std::size_t n_{n};
        std::vector<std::vector<double>> results_;
        std::vector<std::vector<double>> samples_;
        std::vector<ConstraintCoeffs<n>*> cons_ptr_;
        std::unique_ptr<BaseOptimizer<n>> opt_ptr_{nullptr};
};


template<std::size_t n>
template<typename T>
BaseSampler<n>::BaseSampler(const T &lb, const T &ub)
{
    this->setBounds(lb, ub);
};

template<std::size_t n>
void BaseSampler<n>::setBounds(const Vector &lb, const Vector &ub)
{
    lb_=lb;
    ub_=ub;
};


template<std::size_t n>
void BaseSampler<n>::setBounds(const std::vector<double> &lb, const std::vector<double> &ub)
{
    assert (ub.size() == lb.size());
    assert (ub.size() == lb_.size());
    for (std::size_t i = 0; i < lb.size(); ++i)
    {   
        lb_(i) = lb[i];
        ub_(i) = ub[i];
    };
};

template<std::size_t n>
bool BaseSampler<n>::hasOptimizer(){
    return (opt_ptr_ ==  nullptr) ? false : true;
}

template<std::size_t n>
std::vector<double> BaseSampler<n>::optimize(std::vector<double> &seed){
    assert (hasOptimizer());
    return this->opt_ptr_->optimize(seed);
}

template<std::size_t n>
std::pair<Matrix<double,n,1>,Matrix<double,n,1>> BaseSampler<n>::getBounds()
{  
    return std::make_pair(lb_,ub_);
}

template<std::size_t n>
std::vector<std::vector<double>> BaseSampler<n>::results(){
    return results_;
}


template<std::size_t n>
std::vector<std::vector<double>> BaseSampler<n>::samples(){
    return samples_;
}

template<std::size_t n>
void BaseSampler<n>::saveSamples(std::string name){
    std::string new_name;
    new_name = name +"_seeds";
    utils::writeVec2File(samples_,new_name);
}


template<std::size_t n>
void BaseSampler<n>::saveResults(std::string name){
    std::string new_name;
    new_name = name +"_samples";
    if (!results_.empty()){
        utils::writeVec2File(results_,new_name);
    }
}

template<std::size_t n>
void BaseSampler<n>::reset(){
    results_.clear();
    samples_.clear();
}

template<std::size_t n>
void BaseSampler<n>::addConstraints(std::vector<ConstraintCoeffs<n>> & cons){
    for (int i=0; i<cons.size(); i++){
        cons_ptr_.push_back(&cons[i]);
    } 
}

template<std::size_t n>
bool BaseSampler<n>::checkFeasible(const std::vector<double> &x){
    if (cons_ptr_.empty()){
        return true;
    } else {
        return isFeasibleM<n>(x, cons_ptr_);
    }
}

template<std::size_t n>
bool BaseSampler<n>::checkFeasible(const Vector &x){
    std::vector<double> x_vec(n);
    utils::copyEig2Vec(x, x_vec);
    return checkFeasible(x_vec);
}
/*
template<std::size_t n>
void BaseSampler<n>::setOptimizer(BaseOptimizer<n> *opt){
    opt_ptr_ = opt;
    if (opt != nullptr) {
        opt_ptr_->addConstraints(cons_ptr_);
    }
}*/

template<std::size_t n>
void BaseSampler<n>::setOptimizer(std::string method, std::vector<double> &lb, std::vector<double> &ub){
    if (method == "biased"){
        opt_ptr_ = std::make_unique<BiasedOptimizer<n>>(lb, ub);
        opt_ptr_->addConstraints(cons_ptr_);
    }
}


template<std::size_t n> class UniformSampler: public BaseSampler<n>
{
    typedef Matrix<double, n, 1> Vector;
    enum {NeedsToAlign = (sizeof(Vector)%16)==0};

    public:

        UniformSampler(){};
        template <typename T> UniformSampler(const T &lb, const T &ub);
        template <typename T> void setBounds(const T &lb, const T &ub);
        Vector sample();
        Vector sample(const Vector &range, const Vector &lb);
        virtual void run(int n_iter);
        virtual void run(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub);
        void run(const int n_iter, Vector &lb, Vector &ub);

    private:  

        Vector x_temp_;
        Vector range_;
        Vector ones_;
};

template<std::size_t n> 
template<typename T>
UniformSampler<n>::UniformSampler(const T &lb, const T &ub): BaseSampler<n>(lb, ub) 
{
    this->setBounds(lb, ub);
}

template<std::size_t n>
Matrix<double, n, 1> UniformSampler<n>::sample(const Vector& range, const Vector& lb)
{
    this->x_.setRandom();
    x_temp_ = ((this->x_ + ones_) * 0.5).cwiseProduct(range)  +lb; 
    return x_temp_;
}

template<std::size_t n>
Matrix<double, n, 1> UniformSampler<n>::sample()
{
    return this->sample(range_, this->lb_);
}


template<std::size_t n>
template<typename T>
void UniformSampler<n>::setBounds(const T &lb, const T &ub)
{
    BaseSampler<n>::setBounds(lb, ub);
    range_ = this->ub_ - this->lb_;
    ones_.setOnes();
}

template<std::size_t n>
void UniformSampler<n>::run(int n_iter, Vector &lb, Vector &ub){
    std::vector<double> sample_vec(n); 
    std::vector<double> result(n);
    for (int i=0; i<n_iter; i++){
        utils::copyEig2Vec(sample(), sample_vec);
        this->samples_.push_back(sample_vec);
        if (this->opt_ptr_ == nullptr) {
           //if (isFeasibleM<n>(sample_vec, this->cons_ptr_) && boundsCheckVec<n>(sample_vec,lb,ub)) {
           if (this->checkFeasible(sample_vec) && boundsCheckVec<n>(sample_vec,lb,ub)) {
              this->results_.push_back(sample_vec); 
           }
        } else {
           if (boundsCheckVec<n>(sample_vec,lb,ub)){
              result = this->opt_ptr_->optimize(sample_vec);
              this->results_.push_back(result);
           }
        }
    }
}

template<std::size_t n>
void UniformSampler<n>::run(int n_iter){
    run(n_iter, this->lb_, this->ub_);
}

template<std::size_t n>
void UniformSampler<n>::run(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub){
    Matrix<double,n,1> seed_eig(seed.data());
    Matrix<double,n,1> lb_old, ub_old, lb_new(lb.data()), ub_new(ub.data());
    lb_old = this->lb_;
    ub_old = this->ub_;
    this->setBounds(seed_eig + lb_new , seed_eig + ub_new);
    run(n_iter, lb_old, ub_old);
    this->setBounds(lb_old, ub_old);
}

template<std::size_t n> class MetropolisHastings: public BaseSampler<n>
{
    typedef Matrix<double, n, 1> Vector;
    enum {NeedsToAlign = (sizeof(Vector)%16)==0};

    public:

        MetropolisHastings(){};
        template <typename T> MetropolisHastings(const T &lb, const T &ub);//: BaseSampler<n>(lb,ub){};
        //void setP(double (*p)(const Vector&));
        void setP(std::function<double(const Vector&)> &p);
        //void setQ(double (*q)(const Vector&, const Vector&));
        void setA(double (*A)(const Vector&, const Vector&));
        //void set_Q(Vector (*Q)(const Vector&));
        void setQSampler(std::function<Vector(const Vector&)> &Q);
        void setQ(std::function<double(const Vector&, const Vector&)> &q);
        double getP(const Vector &x);
        double getQ(const Vector &x_star, const Vector &x_i);
        Vector sampleQ(const Vector &x);
        void run(const int n_iter);
        double aDefault(const Vector &x_star, const Vector &x_i);
        void saveResults(const std::string &name);
        void saveSamples(const std::string &name);
        std::vector<std::vector<double>> results();
        /* virtual void run(int n_iter){};
        virtual void run(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub){};
        virtual void runOnTangent(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub){};
        *///void set_Q_2(std::function<Vector(const Vector&)> &Q2);

    private:

        bool use_default_A_{true};
        std::vector<std::vector<double>> samples_; // x_star
        std::vector<std::vector<double>> results_; //x_i
        UniformSampler<n> start_;
        int n_samples;
        //double (*p_)(const Vector&);
        //double (*q_)(const Vector&, const Vector&);
        double (*A_)(const Vector&, const Vector&);
        std::function<double(const Vector&, const Vector&)> q_;
        std::function<double(const Vector&)> p_;
        std::function<Vector(const Vector&)> Q_;
        //Vector& (*Q_)(const Vector&, void*data);
};


// this is not finished
// return results propoerly
template<std::size_t n>
std::vector<std::vector<double>> MetropolisHastings<n>::results()
{
    return results_;
};

template<std::size_t n>
template<typename T>
MetropolisHastings<n>::MetropolisHastings(const T &lb, const T &ub)
{
    start_.setBounds(lb, ub);
};

template<std::size_t n>
void MetropolisHastings<n>::setQSampler(std::function<Vector(const Vector&)> &Q)
{
    Q_ = Q;
};

template<std::size_t n>
void MetropolisHastings<n>::setQ(std::function<double(const Vector&, const Vector&)> &q)
{
    q_ = q;
};

template<std::size_t n>
Matrix<double,n,1> MetropolisHastings<n>::sampleQ(const Vector &x)
{
    return Q_(x);
};

    /*
template<std::size_t n>
void MetropolisHastings<n>::set_Q(Vector (*Q)(const Vector&))
{
    this->Q_ = Q;
};*/

template<std::size_t n>
void MetropolisHastings<n>::setP(std::function<double(const Vector&)> &p)
{
    p_ = p;
};
/*
template<std::size_t n>
void MetropolisHastings<n>::setQ(double (*q)(const Vector&, const Vector&))
{
    this->q_ = q;
};*/

template<std::size_t n>
void MetropolisHastings<n>::setA(double (*A)(const Vector&, const Vector&))
{
    this->A_ = A;
}


template<std::size_t n>
double MetropolisHastings<n>::getP(const Vector &x)
{
    return this->p_(x);
};

template<std::size_t n>
double MetropolisHastings<n>::getQ(const Vector &x_star, const Vector &x_i)
{
    return this->q_(x_star, x_i);
};


template<std::size_t n>
double MetropolisHastings<n>::aDefault(const Vector &x_star, const Vector &x_i)
{
    double p_star, q_star, pq_q, q_i, p_i;
    p_star = this->p_(x_star);
    q_star = this->q_(x_star, x_i);
    p_i = this->p_(x_i);
    q_i = this->q_(x_i, x_star);
    pq_q = (p_star/q_star) / (p_i/q_i);
    return (pq_q < 1) ? pq_q : 1; 
}; 




template<std::size_t n>
void MetropolisHastings<n>::run(int n_iter)
{

    std::default_random_engine generator;
    std::uniform_real_distribution<double> u_sampler(0.0,1.0);
    // initialize starting point
    Vector x_i = start_.sample();
    Vector x_star;
    double A, u;
    if (this->use_default_A_) 
    {
        for (int i=0; i< n_iter; i++)
        {
            u = u_sampler(generator);
            x_star = Q_(x_i);
            A = aDefault(x_star, x_i);
            if (u<A)
            {
                x_i = x_star;
                results_.push_back(utils::copyEig2Vec(x_i));
            }               
            samples_.push_back(utils::copyEig2Vec(x_star));
        }
    }
};

template<std::size_t n>
void MetropolisHastings<n>::saveResults(const std::string &name)
{
    utils::writeVec2File(results_, name);
};

template<std::size_t n>
void MetropolisHastings<n>::saveSamples(const std::string &name)
{
    utils::writeVec2File(samples_, name);
};


template <std::size_t n>
class UniformNeighborhoodSampler: public BaseSampler<n>
{
    typedef Matrix<double, n, 1> Vector;
    enum {NeedsToAlign = (sizeof(Vector)%16)==0};

    public:
        UniformNeighborhoodSampler(){};
        template<typename T> UniformNeighborhoodSampler(const T &lb, const T &ub): BaseSampler<n>(lb, ub){};
        // proposal prob density
        double operator()(const Vector& x_star, const Vector& x_i);
        // sample from proposal dist
        Vector operator()(const Vector& x);
        void setWidths(const Vector &widths); 

    private:
        Vector widths_, lb_i_, ub_i_;
        UniformSampler<n> uni_;
};

// target prob dens
template<std::size_t n>
struct TargetProb
{
    std::vector<ConstraintCoeffs<n>> cons;
    double slack{0};
    ConstraintCoeffs<n> cons_temp;
    double operator()(const Matrix<double,n,1> &x){
        typename std::vector<ConstraintCoeffs<n>>::iterator it;
        for(it = cons.begin(); it != cons.end(); ++it)
        {
            cons_temp = *it;
            if (cons_temp.constype == "linear")
            {
                slack = linearConstraint<n>(x,&cons_temp);
            } else if (cons_temp.constype == "quadratic") {
                slack = quadraticConstraint<n>(x,&cons_temp);
            } else {
                std::cout << "wrong constraint type" << std::endl;
            };
            if (slack > 0) return 0;
        }
        return 1;
    }
};

template<std::size_t n>
void UniformNeighborhoodSampler<n>::setWidths(const Vector &widths)
{
    widths_ = widths;
}


// proposal prob dens
template<std::size_t n>
double UniformNeighborhoodSampler<n>::operator()(const Vector &x_star, const Vector &x_i)
{
    return 1;
};

// sample from prop dist
template<std::size_t n>
Matrix<double,n,1> UniformNeighborhoodSampler<n>::operator()(const Vector &x)
{
    lb_i_ = x - widths_;
    ub_i_ = x + widths_; 
    lb_i_ = lb_i_.cwiseMax(this->lb_);
    ub_i_ = ub_i_.cwiseMin(this->ub_);  
    uni_.setBounds(lb_i_, ub_i_);
    return uni_.sample();
};

template<std::size_t n>
class GridWalk : public MetropolisHastings<n> {

    typedef Matrix<double, n, 1> Vector;
    enum {NeedsToAlign = (sizeof(Vector)%16)==0};

    public:
     
        GridWalk(){};
        GridWalk(const std::vector<double> &lb, const std::vector<double> &ub, const::std::vector<double> &widths);
        double PGridWalk(const Matrix<double,n,1> &x);
        Matrix<double,n,1> qSamplerGridWalk(const Matrix<double,n,1> &x);
        Matrix<double,n,1> qSamplerBallWalk(const Matrix<double,n,1> &x);
        double QGridWalk(const Matrix<double,n,1> &x_star, const Matrix<double,n,1> &x);
        void makeBallWalk(double radius);

    protected:

        Matrix<double,n,1> widths_;
        double r_;
    
};

template<std::size_t n>
GridWalk<n>::GridWalk(const std::vector<double> &lb, const std::vector<double> &ub, const std::vector<double> &widths)  {
    this->setBounds(lb,ub);
    widths_ = Matrix<double,n,1>(widths.data()); 
    std::function<double(const Matrix<double,n,1>&)> p = std::bind(&GridWalk<n>::PGridWalk, this, std::placeholders::_1);
    std::function<Matrix<double,n,1>(const Matrix<double,n,1>&)> q_sample = std::bind(&GridWalk<n>::qSamplerGridWalk, this, std::placeholders::_1);
    std::function<double(const Matrix<double,n,1>&, const Matrix<double,n,1>&)> q = std::bind(&GridWalk<n>::QGridWalk, this, std::placeholders::_1, std::placeholders::_2);
    this->setP(p);
    this->setQSampler(q_sample);
    this->setQ(q);
}

template<std::size_t n>
double GridWalk<n>::PGridWalk(const Matrix<double,n,1> &x){
    return (this->checkFeasible(x) && boundsCheck<n>(x, this->lb_, this->ub_)) ? 1 : 0;
}

template<std::size_t n>
Matrix<double,n,1> GridWalk<n>::qSamplerGridWalk(const Matrix<double,n,1> &x){
    Matrix<double,n,1> lb = x - widths_;
    Matrix<double,n,1> ub = x + widths_;
    UniformSampler<n> qsampler(lb, ub);
    return qsampler.sample();
}

template<std::size_t n>
double GridWalk<n>::QGridWalk(const Matrix<double,n,1> &x_star, const Matrix<double,n,1> &x){
    return 1;
}

template<std::size_t n>
Matrix<double,n,1> GridWalk<n>::qSamplerBallWalk(const Matrix<double,n,1> &x){
    Matrix<double,n,1> candidate;
    double d;
    Matrix<double,n,1> lb = x - widths_;
    Matrix<double,n,1> ub = x + widths_;
    UniformSampler<n> qsampler(lb, ub);
    do {
        candidate = qsampler.sample();
        d = (candidate - x).norm();
        std::cout<< "x: \n"  << x << std::endl;
        std::cout<< "c: \n"  << candidate << std::endl;
        std::cout<< "d: \n"  << d << std::endl;
    } while (d > r_) ;
    return candidate;
}

template<std::size_t n>
void GridWalk<n>::makeBallWalk(double radius) {
    r_ = radius;
    std::function<Matrix<double,n,1>(const Matrix<double,n,1>&)> q = std::bind(&GridWalk<n>::qSamplerBallWalk, this, std::placeholders::_1);
    this->setQSampler(q);
}


#endif
