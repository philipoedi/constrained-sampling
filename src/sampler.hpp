#ifndef SAMPLER_H
#define SAMPLER_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <random>
#include <functional>
#include "utils.hpp"
#include "constraints.hpp"

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

    protected:
        Vector lb_;
        Vector ub_;
        Vector x_; 
        std::size_t n_{n};
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
std::pair<Matrix<double,n,1>,Matrix<double,n,1>> BaseSampler<n>::getBounds()
{  
return std::make_pair(lb_,ub_);
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



template<std::size_t n> class MetropolisHastings: public BaseSampler<n>
{
    typedef Matrix<double, n, 1> Vector;
    enum {NeedsToAlign = (sizeof(Vector)%16)==0};

    public:

        MetropolisHastings();
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
        //void set_Q_2(std::function<Vector(const Vector&)> &Q2);

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
            std::cout << A << std::endl;
            if (u<A)
            {
                x_i = x_star;
                results_.push_back(utils::copyEig2Vec(x_i));
            }               
            samples_.push_back(utils::copyEig2Vec(x_star));
        }
    }
        // append data to some datastructure
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
class NormalNeighborhoodSampler : public BaseSampler<n>
{
    typedef Matrix<double, n, 1> Vector;
    enum {NeedsToAlign = (sizeof(Vector)%16)==0};


    public:
        NormalNeighborhoodSampler;
        void setStd(const double sdev);
        double operator()(const Vector &x_star, const Vector &x_i);
        Vector operator()(const Vector &x);
    private:
        Vector sdevs;
}


/*
template<std::size_t n>
class uniform_neighborhood_sampler : public UniformSampler<n>
{

    typedef Matrix<double, n, 1> Vector;
    enum {NeedsToAlign = (sizeof(Vector)%16)==0};

    public:
        uniform_neighborhood_sampler():UniformSampler<n>(){};
        template<typename T> uniform_neighborhood_sampler(const T lb, const T ub)
        : UniformSampler<n>(lb, ub){};
        static Vector sample(const Vector& x);    
        double get_P(const Vector& x);
        void set_P(double (*P)(const Vector&));

        Vector widths_;
    private:
        Vector lb_nb_; // lower bounds for neighborhood
        Vector ub_nb_; // uupper bounds for neighborhood
        Vector range_i_;
        double (*P_)(const Vector &x);
};

template<std::size_t n>
Matrix<double, n, 1> uniform_neighborhood_sampler<n>::sample(const Vector &x)
{
    Vector lb = x - widths_;
    Vector ub = x + widths_;
    Vector range;
    lb_ = lb.cwiseMax(this->lb_);
    ub = lb_.cwiseMin(this->ub_);
    range = ub - lb;
    return this->sample(range, lb);
};

template<std::size_t n>
void uniform_neighborhood_sampler<n>::set_P(double (*P)(const Vector &x))
{
    this->P_ = P;
}

template<std::size_t n>
double uniform_neighborhood_sampler<n>::get_P(const Vector &x)
{
    return P_(x);
}
*/



// dimension

// bounds
 
// uniform sampler

// x_star

// q

//p
//A
#endif
