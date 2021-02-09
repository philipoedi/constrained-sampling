#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;


template<std::size_t n> class uniform_sampler
{
    typedef Matrix<double, n, 1> Vector;
    enum {NeedsToAlign = (sizeof(Vector)%16)==0};

    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign); 
        uniform_sampler(const Vector& lb, const Vector& ub);
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
uniform_sampler<n>::uniform_sampler(const Vector& lb, const Vector& ub)
{
    lb_ = lb;
    ub_ = ub;
    range_ = ub - lb;
    ones_.setOnes();
}

template<std::size_t n>
Matrix<double, n, 1>& uniform_sampler<n>::sample()
{
    x_.setRandom();
    x_temp_ = ((x_ + ones_) * 0.5).cwiseProduct(range_)  + lb_;
    return x_temp_;
}


