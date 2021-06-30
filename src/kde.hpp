#ifndef KDE_H  
#define KDE_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include "utils.hpp"
#include <string>

using namespace Eigen;


template<std::size_t n, std::size_t d> class Kernel
{
    typedef Matrix<double, d, 1> Vector;
//    typedef Matrix<double, n, d> dataMatrix;

    public:
        
        Kernel();
        Kernel(const Vector& bandwidth);
        void fit(const MatrixXd &data);
        double evaluate(const Vector& x);
        void setBandwidth(const Vector& bandwidth);    
        //void find_constant();
        void find_optimal_bandwidth(std::string method);
        double silverman();
        double scott();
        void addData(const std::vector<double> data, const std::size_t i);
        void resize(std::size_t n_rows, std::size_t n_cols);

    private:

        Vector bandwidth_; // 1/bandwidth
  //      dataMatrix data_;
        MatrixXd data_;
  //      dataMatrix distances_;
        double nh_; // num_samples * product over hi 
        std::size_t n_ = n;
        std::size_t d_ = d;
        double C_;
        int nu_ = 2;
        double R_ = 3./5.;
        double kappa_ = 1./5.;
        std::size_t num_rows{n};
               
};

template<std::size_t n, std::size_t d>
Kernel<n,d>::Kernel(){}


template<std::size_t n, std::size_t d>
Kernel<n,d>::Kernel(const Vector& bandwidth)
{
    this->setBandwidth(bandwidth);  
} 

template<std::size_t n, std::size_t d>
void Kernel<n,d>::fit(const MatrixXd &data)
{
    data_.resize(data.rows(), data.cols());
    data_ = data;
}

template<std::size_t n, std::size_t d>
void Kernel<n,d>::addData(const std::vector<double> data, const std::size_t i)
{
    for (int j=0; j<d; j++)
    {
        data_(i,j) = data[j];
    };
}

template<std::size_t n, std::size_t d>
double Kernel<n,d>::evaluate(const Vector& x)
{
    double prob;
    MatrixXd distances;
    distances.resize(data_.rows(),data_.cols());
    distances = (data_- x.transpose().replicate(n,1)).cwiseProduct(bandwidth_.transpose().replicate(n,1));
    distances = (distances.array().abs() > 1.0).select(1, distances);  // select(1 instead of 0 -> 
    distances = 3./4. *( 1.0 - distances.array().square()); // for all vals > 0 follows that distances = 0 because 1-1
    prob = distances.rowwise().prod().sum()*nh_; 
    return prob;
}

template<std::size_t n, std::size_t d>
void Kernel<n,d>::setBandwidth(const Vector& bandwidth)
{
   bandwidth_ = bandwidth; 
   nh_ = bandwidth.prod()/n_;
}

template<std::size_t n, std::size_t d>
void Kernel<n,d>::resize(std::size_t n_rows, std::size_t n_cols){
    data_.resize(n_rows, n_cols);
}
/*
template<std::size_t n, std::size_t d>
void Kernel<n,d>::find_constant()
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
double Kernel<n,d>::silverman()
{
    return pow(4./(d_+2.),(1./(d_+4)))*pow(n_,-1./(d+4.));
}

template<std::size_t n, std::size_t d>
double Kernel<n,d>::scott()
{
    return pow(n_,-1./(d_+4.));
}

template<std::size_t n, std::size_t d>
void Kernel<n,d>::find_optimal_bandwidth(std::string method)
{
    double C;
    assert (method == "scott" || method == "silverman");
    if (method == "scott") C = this->scott();
    else C = this->silverman();
    bandwidth_ = C*((data_.rowwise() - data_.colwise().mean()).array().square().colwise().sum() / (n_-1.)).sqrt();
    std::cout << "bandwidth:" << bandwidth_.array().pow(-1) << std::endl;
    this->setBandwidth(bandwidth_.array().pow(-1)) ;
    //this->setBandwidth(bandwidth_) ;
}




/*template<std::size_t n, std::size_t d>
void Kernel<n,d>::find_constant()
{
}*/


template<std::size_t n, std::size_t d> class KernelEstimator
{
    
    typedef Matrix<double, d, 1> Vector;
    //typedef Matrix<double, n, d> dataMatrix;
    

    public:
        
        KernelEstimator();
        KernelEstimator(const Vector& bandwidth);
        KernelEstimator(const std::string& bandwidth_est);
        void fit(const MatrixXd& data);
        void fit(const std::vector<std::vector<double>>& data);
        void predict(const std::vector<Vector>& x, std::vector<double>& res);
        void predict(const std::vector<double>& lb, const std::vector<double>& ub, const double step);
        void setBandwidth(const Vector& bandwidth);
        void find_optimal_bandwidth(const std::string bandwidth_est);
        void savePdes(const std::string name);
        void evaluateAndAdd(Vector &point);
        void predictOnGrid(const std::vector<double> &lb, const std::vector<double> &ub, const double step);
        void evaluateOnGrid(const std::vector<std::vector<double>> &space, std::size_t vec_index, std::vector<double> &vec_so_far);
        void predictOnSphere(int n_points, double r);
        void setSphere(double r);
        void resizeDataMatrix(std::size_t n_rows);
    
    private:
        
        Kernel<n,d> k_;
        std::string method_{"grid"};
        std::size_t n_{n};
        std::size_t d_{d};
        std::vector<std::vector<double>> pdes_;
        double r_;
};


template<std::size_t n, std::size_t d>
KernelEstimator<n,d>::KernelEstimator(){}

template<std::size_t n, std::size_t d>
KernelEstimator<n,d>::KernelEstimator(const Vector& bandwidth)
{
   k_.setBandwidth(bandwidth);
}

template<std::size_t n, std::size_t d>
void KernelEstimator<n,d>::setBandwidth(const Vector& bandwidth)
{
   k_.setBandwidth(bandwidth);
}

template<std::size_t n, std::size_t d>
void KernelEstimator<n,d>::fit(const MatrixXd &data){
    k_.fit(data);
}

template<std::size_t n, std::size_t d>
KernelEstimator<n,d>::KernelEstimator(const std::string& bandwidth_est)
{
    k_.find_optimal_bandwidth(bandwidth_est);
}

template<std::size_t n, std::size_t d>
void KernelEstimator<n,d>::find_optimal_bandwidth(const std::string bandwidth_est)
{
    k_.find_optimal_bandwidth(bandwidth_est);
}

template<std::size_t n, std::size_t d>
void KernelEstimator<n,d>::setSphere(double r){
    r_ = r;
    method_ = "sphere";
}

template<std::size_t n, std::size_t d>
void KernelEstimator<n,d>::predict(const std::vector<Vector>& x, std::vector<double>& res)
{
    for(std::size_t i = 0; i < x.size(); ++i){
        res[i] = k_.evaluate(x[i]);
    }
}

/*template<std::size_t n, std::size_t d>
void KernelEstimator<n,d>::predict(const std::vector<double>& lb, const std::vector<double>& ub, const double step)
{
    Vector p;
    double x{lb[0]};
    double y{lb[1]};
    double prob;
    pdes_.clear();
    std::vector<double> data;
    data.resize(d+1);
    while (x <= ub[0])
    {
        p(0) = x;
        while (y <= ub[1])
        {
            p(1) = y;
            std::cout << p << std::endl;
            prob = k_.evaluate(p);
            data[0] = x;
            data[1] = y;
            data[2] = prob; 
            y += step;
            pdes_.push_back(data);
        }
        x += step;
        y = lb[1];
    }
}*/

template<std::size_t n, std::size_t d>
void KernelEstimator<n,d>::predict(const std::vector<double> &lb, const std::vector<double> &ub, const double step) {
    assert (lb.size() == ub.size());
    assert (lb.size() == d);
    if (method_ == "grid"){
        predictOnGrid(lb, ub, step);
    } else if (method_ == "sphere") {
        int n_samples{0};
        for (int i=0; i<lb.size(); i++){
            n_samples +=  round((ub[i] - lb[i]) / step);
        }
        predictOnSphere(n_samples, r_);
    }
}

template<std::size_t n, std::size_t d>
void KernelEstimator<n,d>::evaluateAndAdd(Vector &point){
    double prob;
    std::vector<double> data(d);
    prob = k_.evaluate(point);
    utils::copyEig2Vec(point, data); 
    data.push_back(prob);
    pdes_.push_back(data);
}

template<std::size_t n, std::size_t d>
void KernelEstimator<n,d>::predictOnGrid(const std::vector<double> &lb, const std::vector<double> &ub, const double step){
    assert (lb.size() == ub.size());
    assert (lb.size() == d);
    std::vector<std::vector<double>> space(d);
    int j;
    double p;
    for (int i=0; i<space.size(); i++){
        j = 0;
        p = lb[i];
        while (p <= ub[i]) {
            space[i].push_back(p);
            p += step;
        }
    }
    std::vector<double> vec_so_far(d);
    evaluateOnGrid(space, 0, vec_so_far);
}

template<std::size_t n, std::size_t d>
void KernelEstimator<n,d>::evaluateOnGrid(const std::vector<std::vector<double>> &space, std::size_t vec_index, std::vector<double> &vec_so_far){
    if (vec_index >= space.size()){
        Vector point(vec_so_far.data());
        evaluateAndAdd(point);
        return;
    }
    for (std::size_t i=0; i<space[vec_index].size(); i++){
        vec_so_far[vec_index] = space[vec_index][i];
        evaluateOnGrid(space, vec_index+1, vec_so_far);
    }
}

template<std::size_t n, std::size_t d>
void KernelEstimator<n,d>::predictOnSphere(int n_points, double r){
    int n_count{0}, M_theta, M_phi;
    double a, di, theta, phi, d_theta, d_phi;
    Vector point;
    a = 4*M_PI*r*r / n_points;
    di = sqrt(a);
    M_theta = round(M_PI/di);
    d_theta = M_PI/M_theta;
    d_phi = a/d_theta;
    for (int m=0; m<(M_theta-1); m++){
        theta = M_PI*(m+0.5)/M_theta; 
        M_phi = round(2*M_PI*sin(theta/d_theta)) ;
        for (int k=0; k<(M_phi-1); k++){
            phi = 2*M_PI*k/M_phi;
            point = Vector(utils::spherical2cartesian(theta, phi, r).data());
            evaluateAndAdd(point);
            n_count += 1;
        }
    }
}


template<std::size_t n, std::size_t d>
void KernelEstimator<n,d>::fit(const std::vector<std::vector<double>>& data)
{
    resizeDataMatrix(data.size());
    for (int i=0; i<data.size(); i++)
    {
        k_.addData(data[i], i);
    }
}

template<std::size_t n, std::size_t d>
void KernelEstimator<n,d>::savePdes(const std::string name)
{
    utils::writeVec2File(pdes_, name);
}

template<std::size_t n, std::size_t d>
void KernelEstimator<n,d>::resizeDataMatrix(std::size_t n_rows){
    k_.resize(n_rows,d);
}


#endif



