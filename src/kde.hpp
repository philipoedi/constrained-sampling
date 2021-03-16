#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;


template<std::size_t n, std::size_t d> class kernel
{
    typedef Matrix<double, d, 1> Vector;
    typedef Matrix<double, d, n> dataMatrix;

    public:
        
        kernel();
        kernel(const Vector& bandwidth);
        void fit(const dataMatrix& data);
        double evaluate(const Vector& x);
        void set_bandwidth(const Vector& bandwidth);    

    private:

        Vector bandwidth_; // 1/bandwidth,
        dataMatrix data_;
        dataMatrix distances_;
        double nh_; // num_samples * product over hi 
        std::size_t n_ = n;
        std::size_t d_ = d;
       
};

template<std::size_t n, std::size_t d>
kernel<n,d>::kernel(){}


template<std::size_t n, std::size_t d>
kernel<n,d>::kernel(const Vector& bandwidth)
{
    this->set_bandwidth(bandwidth);  
} 

template<std::size_t n, std::size_t d>
void kernel<n,d>::fit(const dataMatrix& data)
{
    data_ = data;
}

template<std::size_t n, std::size_t d>
double kernel<n,d>::evaluate(const Vector& x)
{
    distances_ = (data_ - x.replicate(1,n)).cwiseProduct(bandwidth_.replicate(1,n));
    distances_ = (distances_.array().abs() > 1.0).select(0, distances_); 
    distances_ = 3./4. *( 1.0 - distances_.array().square());
    return distances_.rowwise().prod().sum()*nh_;
}

template<std::size_t n, std::size_t d>
void kernel<n,d>::set_bandwidth(const Vector& bandwidth)
{
   bandwidth_ = bandwidth; 
   nh_ = bandwidth.prod()/n_;
}


template<std::size_t n, std::size_t d> class kernel_estimator
{
    
    typedef Matrix<double, d, 1> Vector;
    typedef Matrix<double, d, n> dataMatrix;

    public:
        
        kernel_estimator();
        kernel_estimator(const Vector& bandwidth);
        void fit(const dataMatrix& data);
        void predict(const std::vector<Vector>& x, std::vector<double>& res);
        void set_bandwidth(const Vector& bandwidth);

    private:
        
        kernel<n,d> k_;
        std::size_t n_;
        std::size_t d_;

};


template<std::size_t n, std::size_t d>
kernel_estimator<n,d>::kernel_estimator(){}

template<std::size_t n, std::size_t d>
kernel_estimator<n,d>::kernel_estimator(const Vector& bandwidth)
{
   k_.set_bandwidth(bandwidth);
}

template<std::size_t n, std::size_t d>
void kernel_estimator<n,d>::set_bandwidth(const Vector& bandwidth)
{
   k_.set_bandwidth(bandwidth);
}

template<std::size_t n, std::size_t d>
void kernel_estimator<n,d>::fit(const dataMatrix& data){
    k_.fit(data);
}

template<std::size_t n, std::size_t d>
void kernel_estimator<n,d>::predict(const std::vector<Vector>& x, std::vector<double>& res)
{
    for(std::size_t i = 0; i < x.size(); ++i){
        res[i] = k_.evaluate(x[i]);
    }
}






