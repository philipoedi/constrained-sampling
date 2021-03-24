#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;


template<std::size_t n> class uniform_sampler
{
    typedef Matrix<double, n, 1> Vector;
    enum {NeedsToAlign = (sizeof(Vector)%16)==0};

    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign); 
        uniform_sampler();
        uniform_sampler(const Vector& lb, const Vector& ub);
        void set_bounds(const Vector& lb, const Vector& ub);
        void set_bounds(const std::vector<double>& lb, const std::vector<double>& ub);
        std::pair<Vector,Vector> get_bounds();
        Vector& sample();

    private:  
    
        Vector lb_;
        Vector ub_;
        Vector x_temp_;
        Vector range_;
        Vector x_;
        Vector ones_;
};
 

template<std::size_t n> 
uniform_sampler<n>::uniform_sampler()
{
}

template<std::size_t n> 
uniform_sampler<n>::uniform_sampler(const Vector& lb, const Vector& ub)
{
    this->set_bounds(lb, ub);
}

template<std::size_t n>
Matrix<double, n, 1>& uniform_sampler<n>::sample()
{
    x_.setRandom();
    x_temp_ = ((x_ + ones_) * 0.5).cwiseProduct(range_)  + lb_;
    return x_temp_;
}

template<std::size_t n>
void uniform_sampler<n>::set_bounds(const Vector& lb, const Vector& ub)
{
    lb_ = lb;
    ub_ = ub;
    range_ = ub - lb;
    ones_.setOnes();
}

template<std::size_t n>
void uniform_sampler<n>::set_bounds(const std::vector<double>& lb, const std::vector<double>& ub)
{
    assert (ub.size() == lb.size());
    assert (ub.size() == lb_.size());
    for (std::size_t i = 0; i < lb.size(); ++i)
    {   
        lb_(i) = lb[i];
        ub_(i) = ub[i];
    };
    range_ = ub_ - lb_;
    ones_.setOnes();
}

template<std::size_t n>
std::pair<Matrix<double,n,1>,Matrix<double,n,1>> uniform_sampler<n>::get_bounds()
{  
    return std::make_pair(lb_,ub_);
}



