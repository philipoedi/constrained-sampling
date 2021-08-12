#ifndef SAMPLER_H
#define SAMPLER_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <random>
#include <functional>
#include <cassert>
#include <memory>
#include "tangent.hpp"

template<std::size_t n>
class BaseOptimizer;

template<std::size_t n>
class BiasedOptimizer;

#include "utils.hpp"
#include "constraints.hpp"
#include "optimizer.hpp"

using namespace Eigen;


template<std::size_t n, std::size_t m>
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
        void saveNumIterations(std::string name);
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


template<std::size_t n, std::size_t m>
template<typename T>
BaseSampler<n,m>::BaseSampler(const T &lb, const T &ub)
{
    this->setBounds(lb, ub);
};

template<std::size_t n, std::size_t m>
void BaseSampler<n,m>::setBounds(const Vector &lb, const Vector &ub)
{
    lb_=lb;
    ub_=ub;
};


template<std::size_t n, std::size_t m>
void BaseSampler<n,m>::setBounds(const std::vector<double> &lb, const std::vector<double> &ub)
{
    assert (ub.size() == lb.size());
    assert (ub.size() == lb_.size());
    for (std::size_t i = 0; i < lb.size(); ++i)
    {   
        lb_(i) = lb[i];
        ub_(i) = ub[i];
    };
};

template<std::size_t n, std::size_t m>
bool BaseSampler<n,m>::hasOptimizer(){
    return (opt_ptr_ ==  nullptr) ? false : true;
}

template<std::size_t n, std::size_t m>
std::vector<double> BaseSampler<n,m>::optimize(std::vector<double> &seed){
    assert (hasOptimizer());
    return this->opt_ptr_->optimize(seed);
}

template<std::size_t n, std::size_t m>
std::pair<Matrix<double,n,1>,Matrix<double,n,1>> BaseSampler<n,m>::getBounds()
{  
    return std::make_pair(lb_,ub_);
}

template<std::size_t n, std::size_t m>
std::vector<std::vector<double>> BaseSampler<n,m>::results(){
    return results_;
}


template<std::size_t n, std::size_t m>
std::vector<std::vector<double>> BaseSampler<n,m>::samples(){
    return samples_;
}

template<std::size_t n, std::size_t m>
void BaseSampler<n,m>::saveSamples(std::string name){
    std::string new_name;
    new_name = name +"_seeds";
    utils::writeVec2File(samples_,new_name);
}


template<std::size_t n, std::size_t m>
void BaseSampler<n,m>::saveResults(std::string name){
    std::string new_name;
    new_name = name +"_samples";
    if (!results_.empty()){
        utils::writeVec2File(results_,new_name);
    }
}

template<std::size_t n, std::size_t m>
void BaseSampler<n,m>::saveNumIterations(std::string name){
    std::string new_name;
    new_name = name +"_num_iterations";
    std::vector<int> num_its;
    num_its = opt_ptr_->getNumIterations();
    if (!num_its.empty()){
        utils::writeVec2File(num_its, new_name);
    }
}


template<std::size_t n, std::size_t m>
void BaseSampler<n,m>::reset(){
    results_.clear();
    samples_.clear();
    if (hasOptimizer()) {
        opt_ptr_->reset();
    }
}

template<std::size_t n, std::size_t m>
void BaseSampler<n,m>::addConstraints(std::vector<ConstraintCoeffs<n>> & cons){
    for (int i=0; i<cons.size(); i++){
        cons_ptr_.push_back(&cons[i]);
    } 
}

template<std::size_t n, std::size_t m>
bool BaseSampler<n,m>::checkFeasible(const std::vector<double> &x){
    if (cons_ptr_.empty()){
        return true;
    } else {
        return isFeasibleM<n>(x, cons_ptr_);
    }
}

template<std::size_t n, std::size_t m>
bool BaseSampler<n,m>::checkFeasible(const Vector &x){
    std::vector<double> x_vec(n);
    utils::copyEig2Vec(x, x_vec);
    return checkFeasible(x_vec);
}
/*
template<std::size_t n, std::size_t m>
void BaseSampler<n,m>::setOptimizer(BaseOptimizer<n> *opt){
    opt_ptr_ = opt;
    if (opt != nullptr) {
        opt_ptr_->addConstraints(cons_ptr_);
    }
}*/

template<std::size_t n, std::size_t m>
void BaseSampler<n,m>::setOptimizer(std::string method, std::vector<double> &lb, std::vector<double> &ub){
    if (method == "biased"){
        opt_ptr_ = std::make_unique<BiasedOptimizer<n>>(lb, ub);
        opt_ptr_->addConstraints(cons_ptr_);
    }
}


template<std::size_t n, std::size_t m>
class UniformSampler: public BaseSampler<n,m>
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

template<std::size_t n, std::size_t m> 
template<typename T>
UniformSampler<n,m>::UniformSampler(const T &lb, const T &ub): BaseSampler<n,m>(lb, ub) 
{
    this->setBounds(lb, ub);
}

template<std::size_t n, std::size_t m>
Matrix<double, n, 1> UniformSampler<n,m>::sample(const Vector& range, const Vector& lb)
{
    this->x_.setRandom();
    x_temp_ = ((this->x_ + ones_) * 0.5).cwiseProduct(range)  +lb; 
    return x_temp_;
}

template<std::size_t n, std::size_t m>
Matrix<double, n, 1> UniformSampler<n,m>::sample()
{
    return this->sample(range_, this->lb_);
}


template<std::size_t n, std::size_t m>
template<typename T>
void UniformSampler<n,m>::setBounds(const T &lb, const T &ub)
{
    BaseSampler<n,m>::setBounds(lb, ub);
    range_ = this->ub_ - this->lb_;
    ones_.setOnes();
}

template<std::size_t n, std::size_t m>
void UniformSampler<n,m>::run(int n_iter, Vector &lb, Vector &ub){
    std::vector<double> sample_vec(n); 
    std::vector<double> result(n);
    int num_samples{0};
    while (num_samples < n_iter) {
    //for (int i=0; i<n_iter; i++){
        utils::copyEig2Vec(sample(), sample_vec);
        this->samples_.push_back(sample_vec);
        if (this->opt_ptr_ == nullptr) {
           //if (isFeasibleM<n>(sample_vec, this->cons_ptr_) && boundsCheckVec<n>(sample_vec,lb,ub)) {
           if (this->checkFeasible(sample_vec) && boundsCheckVec<n>(sample_vec,lb,ub)) {
              this->results_.push_back(sample_vec); 
              num_samples++;
           }
        } else {
           if (boundsCheckVec<n>(sample_vec,lb,ub)){
              result = this->opt_ptr_->optimize(sample_vec);
              this->results_.push_back(result);
              num_samples++;
           }
        }
    }
}

template<std::size_t n, std::size_t m>
void UniformSampler<n,m>::run(int n_iter){
    run(n_iter, this->lb_, this->ub_);
}

template<std::size_t n, std::size_t m>
void UniformSampler<n,m>::run(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub){
    Matrix<double,n,1> seed_eig(seed.data());
    Matrix<double,n,1> lb_old, ub_old, lb_new(lb.data()), ub_new(ub.data());
    lb_old = this->lb_;
    ub_old = this->ub_;
    this->setBounds(seed_eig + lb_new , seed_eig + ub_new);
    run(n_iter, lb_old, ub_old);
    this->setBounds(lb_old, ub_old);
}

class SphereSampler : public UniformSampler<2,1>
{
    public:
        template <typename T> SphereSampler(const T &lb, const T &ub);
        virtual void run (int n_iter);
        void setRadius(double r);
        void setCenter(std::vector<double> x0);

    private:
        double r_{1};
        std::vector<double> x0_{0,0};
};

template<typename T>
SphereSampler::SphereSampler(const T &lb, const T &ub): UniformSampler<2,1>(lb,ub){
}

void SphereSampler::setRadius(double r) {
    r_ = r;
}

void SphereSampler::setCenter(std::vector<double> x0){
    x0_ = x0;
}

void SphereSampler::run(int n_iter) {
    std::vector<double> sample_xyz(3,0);
    std::vector<double> sample_polar(2,0);
    double z;
    double phi;
    Map<Matrix<double,2,1>> sample_polar_eig(sample_polar.data());
    for (int i=0; i<n_iter; i++){
        sample_polar_eig = this->sample();
        z = sample_polar_eig(1);
        phi = sample_polar_eig(0);
        sample_xyz = utils::spherical2cartesianSampled(z, r_, phi);
        sample_xyz[0] += x0_[0];
        sample_xyz[1] += x0_[1];
        sample_xyz[2] += x0_[2];
        this->samples_.push_back(sample_xyz);
        this->results_.push_back(sample_xyz);
    }   
}

class CircleSampler : public UniformSampler<1,1> {

    public:
        template<typename T> CircleSampler(const T &lb, const T &ub);
        virtual void run(int n_iter);
        void setRadius(double r);
        void setCenter(std::vector<double> x0);

    private:
        double r_{1};
        std::vector<double> x0_{0,0};
};


template<typename T>
CircleSampler::CircleSampler(const T &lb, const T &ub) : UniformSampler<1,1>(lb,ub){
}

void CircleSampler::setRadius(double r){
    r_ = r;
}

void CircleSampler::setCenter(std::vector<double> x0){
    x0_ = x0;
}

void CircleSampler::run(int n_iter){
    double phi;
    Matrix<double,1,1> sample_eig;
    std::vector<double> sample_vec(2,0);
    for (int i=0; i<n_iter; i++){
        sample_eig = this->sample();
        phi = sample_eig(0);
        sample_vec = utils::polar2cartesian(phi, r_, x0_);
        this->samples_.push_back(sample_vec);
        this->results_.push_back(sample_vec);
    }
}




class LineSampler : public UniformSampler<2,1> {
    
    public:
        template<typename T> LineSampler(const T &lb, const T &ub);
        void getMap();
        virtual void run(int n_iter);

    public:
        std::vector<double> b1_;
        std::vector<double> b2_;
        double a1_,a2_;
};


template<typename T>
LineSampler::LineSampler(const T &lb, const T &ub) : UniformSampler<2,1>(lb,ub){
}


void LineSampler::getMap(){
    std::vector<std::vector<double>> bounds;
    std::vector<double> coords_on_bound(2,0);
    // iterate over lb
    ConstraintCoeffs<2> *  con = this->cons_ptr_[0];
    double c = con->cons;
    double  x1, x2;
    a1_ = con->coeffs[0];
    a2_ = con->coeffs[1];

    x1 = this->lb_(0);
    x2 = (-c - (a1_*x1)) / a2_;
    coords_on_bound[0] = x1;
    coords_on_bound[1] = x2;
    if (boundsCheckVec<2>(coords_on_bound, this->lb_, this->ub_)){
        bounds.push_back(coords_on_bound);
    }
      
    x1 = this->ub_(0);
    x2 = (-c - (a1_*x1)) / a2_;
    coords_on_bound[0] = x1;
    coords_on_bound[1] = x2;
    if (boundsCheckVec<2>(coords_on_bound, this->lb_, this->ub_)){
        bounds.push_back(coords_on_bound);
    }

    x2 = this->lb_(1);
    x1 = (-c - (a2_*x2)) / a1_;
    coords_on_bound[0] = x1;
    coords_on_bound[1] = x2;
    if (boundsCheckVec<2>(coords_on_bound, this->lb_, this->ub_)){
        bounds.push_back(coords_on_bound);
    }
 

    x2 = this->ub_(1);
    x1 = (-c - (a2_*x2)) / a1_;
    coords_on_bound[0] = x1;
    coords_on_bound[1] = x2;
    if (boundsCheckVec<2>(coords_on_bound, this->lb_, this->ub_)){
        bounds.push_back(coords_on_bound);
    }

    b1_ = bounds[0];
    bool new_bound_found = false;
    int j = 1;
    double sub_space_bound;
    while (!new_bound_found) {
        b2_ = bounds[j];
        if (b2_[0] == b1_[0] && b2_[1] == b2_[1]){
            ++j; 
        } else {
            a1_ = b2_[0] - b1_[0];
            a2_ = b2_[1] - b1_[1];
            sub_space_bound = (b2_[0] - b1_[0]) / a1_;
            new_bound_found = true;
        }
    }
    


    std::vector<double> lb_new(2,0);
    std::vector<double> ub_new(2,0);
    
    if (sub_space_bound < 0) {
        lb_new[0] = sub_space_bound;
        ub_new[0] = 0;
    } else {
        lb_new[0] = 0;
        ub_new[0] = sub_space_bound;
    }
    this->setBounds(lb_new, ub_new);
}


void LineSampler::run(int n_iter){
    getMap();
    std::cout << "a1: " << a1_ << " a2: " << a2_ << std::endl;
    std::cout << "lb: " << this->lb_ << std::endl; 
    std::cout << "ub: " << this->ub_ << std::endl; 
    std::cout << "b1: " << b1_[0] << " " <<b1_[1]  << std::endl; 
    std::cout << "b2: " << b2_[0] << " " <<b2_[1]  << std::endl; 

    double t;
    std::vector<double> sample(2);
    this->cons_ptr_[0];
    for (int i=0; i<n_iter; i++){
        t = this->sample()(0);
        std::cout << "t " << t << std::endl;
        sample[0] = b1_[0] + t * a1_;
        sample[1] = b1_[1] + t * a2_;
        std::cout << sample[0] << " " << sample[1] << std::endl;
    }
    
}

template<std::size_t n, std::size_t m> 
class MetropolisHastings: public BaseSampler<n,m>
{
    typedef Matrix<double, n, 1> Vector;
    enum {NeedsToAlign = (sizeof(Vector)%16)==0};

    public:

        MetropolisHastings(){};
        template <typename T> MetropolisHastings(const T &lb, const T &ub);//: BaseSampler<n,m>(lb,ub){};
        //void setP(double (*p)(const Vector&));
        template<typename T> void setBounds(const T &lb, const T &ub);
        void setP(std::function<double(const Vector&)> &p);
        //void setQ(double (*q)(const Vector&, const Vector&));
        void setA(double (*A)(const Vector&, const Vector&));
        //void set_Q(Vector (*Q)(const Vector&));
        void setQSampler(std::function<Vector(const Vector&)> &Q);
        void setQ(std::function<double(const Vector&, const Vector&)> &q);
        double getP(const Vector &x);
        double getQ(const Vector &x_star, const Vector &x_i);
        Vector sampleQ(const Vector &x);
        void runOnTangent(const int n_iter);
        virtual void run(const int n_iter);
        virtual void run(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub);
        virtual void runOnTangent(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub);
        double aDefault(const Vector &x_star, const Vector &x_i);
        //void saveResults(const std::string &name);
        //void saveSamples(const std::string &name);
        //std::vector<std::vector<double>> results();
        void setProjectOnEachStep(bool proj);
        //virtual void run(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub){};
        //virtual void runOnTangent(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub){};
        //void set_Q_2(std::function<Vector(const Vector&)> &Q2);

    protected:

        bool use_default_A_{true};
        bool project_on_each_step_{true};
 //       std::vector<std::vector<double>> samples_; // x_star
   //     std::vector<std::vector<double>> results_; //x_i
        UniformSampler<n,m> start_;
        int n_samples;
        TangentSpace<n,m> tang_;
        //double (*p_)(const Vector&);
        //double (*q_)(const Vector&, const Vector&);
        double (*A_)(const Vector&, const Vector&);
        std::function<double(const Vector&, const Vector&)> q_;
        std::function<double(const Vector&)> p_;
        std::function<Vector(const Vector&)> Q_;
        //Vector& (*Q_)(const Vector&, void*data);
};

/*
// this is not finished
// return results propoerly
template<std::size_t n, std::size_t m>
std::vector<std::vector<double>> MetropolisHastings<n,m>::results()
{
    return this->results_;
};
*/
template<std::size_t n, std::size_t m>
template<typename T>
MetropolisHastings<n,m>::MetropolisHastings(const T &lb, const T &ub)
{
    setBounds(lb,ub);
    //tart_.setBounds(lb, ub);
};

template<std::size_t n, std::size_t m>
template<typename T>
void MetropolisHastings<n,m>::setBounds(const T &lb, const T &ub){
    BaseSampler<n,m>::setBounds(lb,ub);
    start_.setBounds(lb,ub);
}

template<std::size_t n, std::size_t m>
void MetropolisHastings<n,m>::setQSampler(std::function<Vector(const Vector&)> &Q)
{
    Q_ = Q;
};

template<std::size_t n, std::size_t m>
void MetropolisHastings<n,m>::setQ(std::function<double(const Vector&, const Vector&)> &q)
{
    q_ = q;
};

template<std::size_t n, std::size_t m>
Matrix<double,n,1> MetropolisHastings<n,m>::sampleQ(const Vector &x)
{
    return Q_(x);
};

template<std::size_t n, std::size_t m>
void MetropolisHastings<n,m>::setProjectOnEachStep(bool proj) {
    project_on_each_step_ = proj;
}

    /*
template<std::size_t n, std::size_t m>
void MetropolisHastings<n,m>::set_Q(Vector (*Q)(const Vector&))
{
    this->Q_ = Q;
};*/

template<std::size_t n, std::size_t m>
void MetropolisHastings<n,m>::setP(std::function<double(const Vector&)> &p)
{
    p_ = p;
};
/*
template<std::size_t n, std::size_t m>
void MetropolisHastings<n,m>::setQ(double (*q)(const Vector&, const Vector&))
{
    this->q_ = q;
};*/

template<std::size_t n, std::size_t m>
void MetropolisHastings<n,m>::setA(double (*A)(const Vector&, const Vector&))
{
    this->A_ = A;
}


template<std::size_t n, std::size_t m>
double MetropolisHastings<n,m>::getP(const Vector &x)
{
    return this->p_(x);
};

template<std::size_t n, std::size_t m>
double MetropolisHastings<n,m>::getQ(const Vector &x_star, const Vector &x_i)
{
    return this->q_(x_star, x_i);
};


template<std::size_t n, std::size_t m>
double MetropolisHastings<n,m>::aDefault(const Vector &x_star, const Vector &x_i)
{
    double p_star, q_star, pq_q, q_i, p_i;
    p_star = this->p_(x_star);
    q_star = this->q_(x_star, x_i);
    p_i = this->p_(x_i);
    q_i = this->q_(x_i, x_star);
    pq_q = (p_star/q_star) / (p_i/q_i);
    return (pq_q < 1) ? pq_q : 1; 
}; 



/*
template<std::size_t n, std::size_t m>
void MetropolisHastings<n,m>::run(int n_iter)
{

    //std::default_random_engine generator;
    //std::uniform_real_distribution<double> u_sampler(0.0,1.0);
    // initialize starting point
    std::vector<double> x_star_vec(n,0), x_i_vec(n,0);
    Map<Vector> x_star(x_star_vec.data(),n);
    Map<Vector> x_i(x_i_vec.data(),n);
    x_i = start_.sample();
    //double A, u;
    if (this->use_default_A_) 
    {
        for (int i=0; i< n_iter; i++)
        {
            x_star = Q_(x_i);
            if (this->hasOptimizer()) {
                samples_.push_back(x_star_vec);
                if (boundsCheck<n>(x_star, this->lb_, this->ub_)) {
                    if (!this->checkFeasible(x_star)) {
                        x_i_vec = this->optimize(x_star_vec);
                        results_.push_back(x_i_vec);
                    }
                }
            } else {
                u = u_sampler(generator);
                A = aDefault(x_star, x_i);
                if (u<A)
                {
                    x_i = x_star;
                    results_.push_back(utils::copyEig2Vec(x_i));
                }               
                samples_.push_back(utils::copyEig2Vec(x_star));
            }
        }
    }
};
*/

template<std::size_t n, std::size_t m>
void MetropolisHastings<n,m>::run(int n_iter) {
    std::vector<double> x_star_vec(n,0), x_i_vec(n,0);
    Map<Vector> x_i(x_i_vec.data(),n);
    Map<Vector> x_star(x_star_vec.data(),n);
    for (int i=0; i<3; i++){
        std::cout << start_.sample() << std::endl;
    }
    x_i = start_.sample();
    run(n_iter, x_i_vec, x_i_vec, x_i_vec);
}

template<std::size_t n, std::size_t m>
void MetropolisHastings<n,m>::run(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub) {
    std::vector<double> x_i_vec = seed;
    std::vector<double> x_star_vec(n,0);
    Map<Vector> x_i(x_i_vec.data(),n);
    Map<Vector> x_star(x_star_vec.data(),n);
    for (int i=0; i<n_iter; i++){
        x_star = Q_(x_i);
        this->samples_.push_back(x_star_vec);
        if (boundsCheck<n>(x_star, this->lb_, this->ub_)){
            if (!this->checkFeasible(x_star)) {
                if (this->hasOptimizer()){
                    x_i_vec = this->optimize(x_star_vec);                  
                    this->results_.push_back(x_i_vec);
                }
             } else {
                this->results_.push_back(x_i_vec);
                x_i = x_star;
            }
        }
    }
}


template<std::size_t n, std::size_t m>
void MetropolisHastings<n,m>::runOnTangent(int n_iter){
    std::vector<double> x_i_vec(n,0);
    Map<Vector> x_i(x_i_vec.data(),n);
    x_i = start_.sample();
    x_i_vec = this->optimize(x_i_vec);
    runOnTangent(n_iter, x_i_vec, x_i_vec, x_i_vec);
}

template<std::size_t n, std::size_t m>
void MetropolisHastings<n,m>::runOnTangent(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub) {
    std::vector<double> x_i_vec = seed;
    std::vector<double> x_star_vec(n,0);
    Map<Vector> x_i(x_i_vec.data(), n);
    Map<Vector> x_star(x_star_vec.data(), n);
    if (project_on_each_step_){
        for (int i=0; i<n_iter; i++){
            this->tang_ = tangentSpaceFromConstraints<n,m>(this->cons_ptr_, x_i_vec, 1e-6);
            x_star = Q_(x_i);
            this->samples_.push_back(x_star_vec);
            if (boundsCheck<n>(x_star, this->lb_, this->ub_)){
                x_i_vec = this->optimize(x_star_vec);
                this->results_.push_back(x_i_vec);
           } 
        }
    }
}

/*
template<std::size_t n, std::size_t m>
void MetropolisHastings<n,m>::saveResults(const std::string &name)
{
    utils::writeVec2File(results_, name);
};

template<std::size_t n, std::size_t m>
void MetropolisHastings<n,m>::saveSamples(const std::string &name)
{
    utils::writeVec2File(samples_, name);
};
*/

template <std::size_t n, std::size_t m>
class UniformNeighborhoodSampler: public BaseSampler<n,m>
{
    typedef Matrix<double, n, 1> Vector;
    enum {NeedsToAlign = (sizeof(Vector)%16)==0};

    public:
        UniformNeighborhoodSampler(){};
        template<typename T> UniformNeighborhoodSampler(const T &lb, const T &ub): BaseSampler<n,m>(lb, ub){};
        // proposal prob density
        double operator()(const Vector& x_star, const Vector& x_i);
        // sample from proposal dist
        Vector operator()(const Vector& x);
        void setWidths(const Vector &widths); 

    private:
        Vector widths_, lb_i_, ub_i_;
        UniformSampler<n,m> uni_;
};

// target prob dens
template<std::size_t n, std::size_t m>
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

template<std::size_t n, std::size_t m>
void UniformNeighborhoodSampler<n,m>::setWidths(const Vector &widths)
{
    widths_ = widths;
}


// proposal prob dens
template<std::size_t n, std::size_t m>
double UniformNeighborhoodSampler<n,m>::operator()(const Vector &x_star, const Vector &x_i)
{
    return 1;
};

// sample from prop dist
template<std::size_t n, std::size_t m>
Matrix<double,n,1> UniformNeighborhoodSampler<n,m>::operator()(const Vector &x)
{
    lb_i_ = x - widths_;
    ub_i_ = x + widths_; 
    lb_i_ = lb_i_.cwiseMax(this->lb_);
    ub_i_ = ub_i_.cwiseMin(this->ub_);  
    uni_.setBounds(lb_i_, ub_i_);
    return uni_.sample();
};

template<std::size_t n, std::size_t m>
class GridWalk : public MetropolisHastings<n,m> {

    typedef Matrix<double, n, 1> Vector;
    enum {NeedsToAlign = (sizeof(Vector)%16)==0};

    public:
     
        GridWalk(){};
        GridWalk(const std::vector<double> &lb, const std::vector<double> &ub, const::std::vector<double> &widths);
        GridWalk(const std::vector<double> &lb, const std::vector<double> &ub, const::std::vector<double> &widths, bool use_tangent);
        GridWalk(const std::vector<double> &lb, const std::vector<double> &ub, const::std::vector<double> &widths, bool use_tangent, double r_ball_walk);
        template<typename T> void setBounds(const T &lb, const T &ub);
        double PGridWalk(const Matrix<double,n,1> &x);
        Matrix<double,n,1> qSamplerGridWalk(const Matrix<double,n,1> &x);
        Matrix<double,n,1> qSamplerGridWalkOnTangent(const Matrix<double,n,1> &x);
        Matrix<double,n,1> qSamplerBallWalk(const Matrix<double,n,1> &x);
        Matrix<double,n,1> qSamplerBallWalkOnTangent(const Matrix<double,n,1> &x);
        double QGridWalk(const Matrix<double,n,1> &x_star, const Matrix<double,n,1> &x);
        void makeBallWalk(double radius);
        void setRunOnTangent(bool run);
        virtual void runOnTangent(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub); 

    protected:

        Matrix<double,n,1> widths_;
        double r_;
        bool run_on_tangent_{false}, ball_walk_{false}; 
};

template<std::size_t n, std::size_t m>
GridWalk<n,m>::GridWalk(const std::vector<double> &lb, const std::vector<double> &ub, const std::vector<double> &widths)  {
    //MetropolisHastings(lb, ub);
    setBounds(lb,ub);
    widths_ = Matrix<double,n,1>(widths.data()); 
    std::function<double(const Matrix<double,n,1>&)> p = std::bind(&GridWalk<n,m>::PGridWalk, this, std::placeholders::_1);
    std::function<Matrix<double,n,1>(const Matrix<double,n,1>&)> q_sample = std::bind(&GridWalk<n,m>::qSamplerGridWalk, this, std::placeholders::_1);
    std::function<double(const Matrix<double,n,1>&, const Matrix<double,n,1>&)> q = std::bind(&GridWalk<n,m>::QGridWalk, this, std::placeholders::_1, std::placeholders::_2);
    this->setP(p);
    this->setQSampler(q_sample);
    this->setQ(q);
}

template<std::size_t n, std::size_t m>
template<typename T>
void GridWalk<n,m>::setBounds(const T &lb, const T &ub){
    MetropolisHastings<n,m>::setBounds(lb, ub);
}


template<std::size_t n, std::size_t m>
GridWalk<n,m>::GridWalk(const std::vector<double> &lb, const std::vector<double> &ub, const::std::vector<double> &widths, bool use_tangent):
    GridWalk<n,m>(lb, ub, widths) {
    setRunOnTangent(use_tangent);
}

template<std::size_t n, std::size_t m>
GridWalk<n,m>::GridWalk(const std::vector<double> &lb, const std::vector<double> &ub, const::std::vector<double> &widths, bool use_tangent, double r_ball_walk)  :  
    GridWalk<n,m>(lb, ub, widths, use_tangent) {
    makeBallWalk(r_ball_walk);
}

template<std::size_t n, std::size_t m>
double GridWalk<n,m>::PGridWalk(const Matrix<double,n,1> &x){
    return (this->checkFeasible(x) && boundsCheck<n>(x, this->lb_, this->ub_)) ? 1 : 0;
}

template<std::size_t n, std::size_t m>
Matrix<double,n,1> GridWalk<n,m>::qSamplerGridWalk(const Matrix<double,n,1> &x){
    Matrix<double,n,1> lb = x - widths_;
    Matrix<double,n,1> ub = x + widths_;
    UniformSampler<n,m> qsampler(lb, ub);
    return qsampler.sample();
}

template<std::size_t n, std::size_t m>
Matrix<double,n,1> GridWalk<n,m>::qSamplerGridWalkOnTangent(const Matrix<double,n,1> &x){
    this->tang_.setSamplerBounds((widths_/-2).head(n-m), (widths_/2).head(n-m));
    return this->tang_.sampleInAmbient();
}

template<std::size_t n, std::size_t m>
Matrix<double,n,1> GridWalk<n,m>::qSamplerBallWalkOnTangent(const Matrix<double,n,1> &x){
    double d;
    Matrix<double,n,1> candidate;
    this->tang_.setSamplerBounds((widths_/-2).head(n-m), (widths_/2).head(n-m));
    do {
        candidate = this->tang_.sampleInAmbient();
        d = (candidate - x).norm();
    } while (d < r_);
    return candidate;
}

template<std::size_t n, std::size_t m>
double GridWalk<n,m>::QGridWalk(const Matrix<double,n,1> &x_star, const Matrix<double,n,1> &x){
    return 1;
}

template<std::size_t n, std::size_t m>
Matrix<double,n,1> GridWalk<n,m>::qSamplerBallWalk(const Matrix<double,n,1> &x){
    Matrix<double,n,1> candidate;
    double d;
    Matrix<double,n,1> lb = x - widths_;
    Matrix<double,n,1> ub = x + widths_;
    UniformSampler<n,m> qsampler(lb, ub);
    do {
        candidate = qsampler.sample();
        d = (candidate - x).norm();
    } while (d > r_) ;
    return candidate;
}

template<std::size_t n, std::size_t m>
void GridWalk<n,m>::makeBallWalk(double radius) {
    ball_walk_ = true;
    r_ = radius;
    std::function<Matrix<double,n,1>(const Matrix<double,n,1>&)> q;
    if (run_on_tangent_) {
        q = std::bind(&GridWalk<n,m>::qSamplerBallWalkOnTangent, this, std::placeholders::_1);
    } else {
        q = std::bind(&GridWalk<n,m>::qSamplerBallWalk, this, std::placeholders::_1);
    }
    this->setQSampler(q);
}

template<std::size_t n, std::size_t m>
void GridWalk<n,m>::setRunOnTangent(bool run){
    run_on_tangent_ = run;
    std::function<Matrix<double,n,1>(const Matrix<double,n,1>&)> q;
    if (run) {
        if (ball_walk_) {
            q = std::bind(&GridWalk<n,m>::qSamplerBallWalkOnTangent, this, std::placeholders::_1);
        } else {
            q = std::bind(&GridWalk<n,m>::qSamplerGridWalkOnTangent, this, std::placeholders::_1);
        }
    } else {
         if (ball_walk_) {
            q = std::bind(&GridWalk<n,m>::qSamplerBallWalk, this, std::placeholders::_1);
        } else {
            q = std::bind(&GridWalk<n,m>::qSamplerGridWalk, this, std::placeholders::_1);
        }
    }
    this->setQSampler(q);
}

template<std::size_t n, std::size_t m>
void GridWalk<n,m>::runOnTangent(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub) {
    MetropolisHastings<n,m>::runOnTangent(n_iter, seed, lb, ub);
}



#endif
