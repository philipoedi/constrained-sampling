#ifndef KDE_H  
#define KDE_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include "utils.hpp"
#include <string>

using namespace Eigen;


template<std::size_t n, std::size_t d> class kernel
{
    typedef Matrix<double, d, 1> Vector;
    typedef Matrix<double, n, d> dataMatrix;

    public:
        
        kernel();
        kernel(const Vector& bandwidth);
        void fit(const dataMatrix& data);
        double evaluate(const Vector& x);
        void set_bandwidth(const Vector& bandwidth);    
        //void find_constant();
        void find_optimal_bandwidth(std::string method);
        double silverman();
        double scott();
        void add_data(const std::vector<double> data, const std::size_t i);
    
    private:

        Vector bandwidth_; // 1/bandwidth
        dataMatrix data_;
        dataMatrix distances_;
        double nh_; // num_samples * product over hi 
        std::size_t n_ = n;
        std::size_t d_ = d;
        double C_;
        int nu_ = 2;
        double R_ = 3./5.;
        double kappa_ = 1./5.;
               
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
void kernel<n,d>::add_data(const std::vector<double> data, const std::size_t i)
{
    for (int j=0; j<d; j++)
    {
        data_(i,j) = data[j];
    };
}

template<std::size_t n, std::size_t d>
double kernel<n,d>::evaluate(const Vector& x)
{
    double prob;
    std::cout << "data_: \n" << data_ << std::endl;
    std::cout << "x: \n" << x << std::endl;
    std::cout << "data_ - x.transpose().repliacete(n,1): \n" << data_ - x.transpose().replicate(n,1) <<std::endl;
    std::cout << "bandwidth transpose replicate: \n" << bandwidth_.transpose().replicate(n,1)<<std::endl;
    distances_ = (data_- x.transpose().replicate(n,1)).cwiseProduct(bandwidth_.transpose().replicate(n,1));
    std::cout << "d-x.t*b" << distances_ << std::endl;
    std::cout << ">1: " << distances_ << std::endl;
    distances_ = (distances_.array().abs() > 1.0).select(1, distances_);  // select(1 instead of 0 -> 
    distances_ = 3./4. *( 1.0 - distances_.array().square()); // for all vals > 0 follows that distances = 0 because 1-1
    std::cout << "distances_ : " << distances_ << std::endl;
    std::cout << "distances_ filter: " << distances_ << std::endl;
    prob = distances_.rowwise().prod().sum()*nh_; 
    std::cout << "dist prod:" << distances_.rowwise().prod() << std::endl; 
    return prob;
}

template<std::size_t n, std::size_t d>
void kernel<n,d>::set_bandwidth(const Vector& bandwidth)
{
   bandwidth_ = bandwidth; 
   nh_ = bandwidth.prod()/n_;
}

/*
template<std::size_t n, std::size_t d>
void kernel<n,d>::find_constant()
{
    double num;
    double den;
    double ff2;
    double ff1;
    double exp;
    num = pow(M_PI,(d_/2.)) * pow(2.,(nu_+d_-1.)) * pow(utils::factorial(nu_),2.)* pow(R_, d_);
    ff2 = utils::factorial(utils::factorial(2.*nu_ - 1.));
    ff1 = utils::factorial(utils::factorial(nu_ - 1.));
    den = nu_* pow(kappa_,2.) * (ff2 + (d_-1.)*pow(ff1,2.));
    exp = 1./(2*nu_+d_);
    std::cout << pow(num/den,(1./(2.*nu_ +d_))) << std::endl;
    std::cout << "hi" << std::endl;
    std::cout <<"num: " << num << std::endl; 
    std::cout <<"den: " << den << std::endl; 
    std::cout <<"exp: " << exp << std::endl;
}*/

template<std::size_t n, std::size_t d>
double kernel<n,d>::silverman()
{
    return pow(4./(d_+2.),(1./(d_+4)))*pow(n_,-1./(d+4.));
}

template<std::size_t n, std::size_t d>
double kernel<n,d>::scott()
{
    return pow(n_,-1./(d_+4.));
}

template<std::size_t n, std::size_t d>
void kernel<n,d>::find_optimal_bandwidth(std::string method)
{
    double C;
    assert (method == "scott" || method == "silverman");
    if (method == "scott") C = this->scott();
    else C = this->silverman();
    bandwidth_ = C*((data_.rowwise() - data_.colwise().mean()).array().square().colwise().sum() / (n_-1.)).sqrt();
    std::cout << "bandwidth:" << bandwidth_.array().pow(-1) << std::endl;
    this->set_bandwidth(bandwidth_.array().pow(-1)) ;
    //this->set_bandwidth(bandwidth_) ;
}




/*template<std::size_t n, std::size_t d>
void kernel<n,d>::find_constant()
{
}*/


template<std::size_t n, std::size_t d> class kernel_estimator
{
    
    typedef Matrix<double, d, 1> Vector;
    typedef Matrix<double, n, d> dataMatrix;

    public:
        
        kernel_estimator();
        kernel_estimator(const Vector& bandwidth);
        kernel_estimator(const std::string& bandwidth_est);
        void fit(const dataMatrix& data);
        void fit(const std::vector<std::vector<double>>& data);
        void predict(const std::vector<Vector>& x, std::vector<double>& res);
        void predict(const std::vector<double>& lb, const std::vector<double>& ub, const double step);
        void set_bandwidth(const Vector& bandwidth);
        void find_optimal_bandwidth(const std::string bandwidth_est);


    private:
        
        kernel<n,d> k_;
        std::size_t n_{n};
        std::size_t d_{d};

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
kernel_estimator<n,d>::kernel_estimator(const std::string& bandwidth_est)
{
    k_.find_optimal_bandwidth(bandwidth_est);
}

template<std::size_t n, std::size_t d>
void kernel_estimator<n,d>::find_optimal_bandwidth(const std::string bandwidth_est)
{
    k_.find_optimal_bandwidth(bandwidth_est);
}


template<std::size_t n, std::size_t d>
void kernel_estimator<n,d>::predict(const std::vector<Vector>& x, std::vector<double>& res)
{
    for(std::size_t i = 0; i < x.size(); ++i){
        res[i] = k_.evaluate(x[i]);
    }
}

template<std::size_t n, std::size_t d>
void kernel_estimator<n,d>::predict(const std::vector<double>& lb, const std::vector<double>& ub, const double step)
{
    Vector p;
    double x{lb[0]};
    double y{lb[1]};
    while (x <= ub[0])
    {
        p(0) = x;
        while (y <= ub[1])
        {
            p(1) = y;
            std::cout << p << std::endl;
            std::cout << "prob: " << k_.evaluate(p)  <<"\n eval end"  << std::endl;
            
            y += step;
        }
        x += step;
        y = lb[1];
    }
}


template<std::size_t n, std::size_t d>
void kernel_estimator<n,d>::fit(const std::vector<std::vector<double>>& data)
{
   for (int i=0; i<data.size(); i++)
    {
        k_.add_data(data[i], i);
    }
}
#endif



