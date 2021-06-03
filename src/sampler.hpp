#ifndef SAMPLER_H
#define SAMPLER_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <random>
#include <functional>

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
        void setA(double (*A)(const double&, const Vector&, const Vector&));
        //void set_Q(Vector (*Q)(const Vector&));
        void setQ(std::function<Vector(const Vector&)> &Q);
        double getP(const Vector &x);
        double getQ(const Vector &x_star, const Vector &x_i);
        Vector sampleQ(const Vector &x);
        void run(const int n_iter);
        double aDefault(const double &pq_i, const Vector &x_star, const Vector &x_i);
        //void set_Q_2(std::function<Vector(const Vector&)> &Q2);

    private:

        bool use_default_A_{true};
        UniformSampler<n> start_;
        int n_samples;
        Vector x_i_;
        Vector x_star_;
        double p_star_;
        double q_star_;
        double p_i_;
        double q_i_;
        double pq_i_;
        //double (*p_)(const Vector&);
        double (*q_)(const Vector&, const Vector&);
        double (*A_)(const double&, const Vector&, const Vector&);
        std::function<double(const Vector&)> p_;
        std::function<Vector(const Vector&)> Q_;
        //Vector& (*Q_)(const Vector&, void*data);
};

template<std::size_t n>
MetropolisHastings<n>::MetropolisHastings(): BaseSampler<n>()
{
};

template<std::size_t n>
template<typename T>
MetropolisHastings<n>::MetropolisHastings(const T &lb, const T &ub)
{
    start_.setBounds(lb, ub);
};

template<std::size_t n>
void MetropolisHastings<n>::setQ(std::function<Vector(const Vector&)> &Q)
{
    Q_ = Q;
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
void MetropolisHastings<n>::setA(double (*A)(const double&, const Vector&, const Vector&))
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
double MetropolisHastings<n>::aDefault(const double &pq_i, const Vector &x_star, const Vector &x_i)
{
    double p_star, q_star, pq_q;
    p_star = this->p_(x_star);
    q_star = this->q_(x_star, x_i);
    pq_q = (p_star/q_star) / (pq_i);
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
            A = aDefault(pq_i, x_star, x_i);
            x_i = (u < A) ? x_star : x_i;
            std::cout << x_i << std::endl;
        }
    }
        // append data to some datastructure
};


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
