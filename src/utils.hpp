#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <algorithm>
#include <random>
#include <chrono>

using namespace Eigen;
/*
string combine_sample_result();

void write_to_file(ofstream file, string results);


void run_sampling();
*/
namespace utils
{
    int factorial(int n);
    
    void copyVec2Vec(const std::vector<double>& src, std::vector<double>& dst);

    void copyEig2Vec(const VectorXd& src, std::vector<double>& dst);
    
    void copyEig2Arr(const VectorXd& src, double *dst);    

    std::vector<double> copyEig2Vec(const VectorXd& src);

    void copyMatvec2Matvec(const std::vector<std::vector<double>>& src, std::vector<std::vector<double>>& dst);

    void writeVec2File(const std::vector<std::vector<double>> &data, const std::string &name);
    
    void writeVec2File(const std::vector<int> &data, const std::string &name);
    void writeVec2File(const std::vector<double> &data, const std::string &name);

    void writeVec2File(const std::vector<std::vector<int>> &data, const std::string &name);
   
//    void write2file_results(const std::vector<std::vector<double>>& seeds, const std::vector<std::vector>>& samples);
    void writeMetadata2File(const std::string &problem_file, const std::string &name);
    
    void appendVec2Vec(const std::vector<std::vector<double>> &src, std::vector<std::vector<double>> &dst);

    std::vector<double> slice(const std::vector<double> &x, int start, int end);

    std::string getDateString();
    std::string getDateTimeString();

    std::vector<double> spherical2cartesian(double theta, double phi, double r);
    std::vector<double> spherical2cartesianSampled(double z, double r,  double phi);
    std::vector<double> polar2cartesian(double phi, double r, std::vector<double> x0);


    void shuffleVector(std::vector<std::vector<double>> &data);
  

    template<std::size_t n>
    double squaredDist(const Matrix<double,n,1> &v1, const Matrix<double,n,1> &v2) 
    { 
       return (v1 - v2).squaredNorm(); 
    };
    
    MatrixXd MatVec2EigMat(const std::vector<std::vector<double>> &data);

    template<std::size_t n>
    Matrix<double,n,1> sampleMean(const std::vector<std::vector<double>> &samples){
        int n_dim, m_dim;
        n_dim = samples.size();
        m_dim = samples[0].size();
        MatrixXd data(n_dim,m_dim);
        data = MatVec2EigMat(samples);
        return data.colwise().mean();
       };
    
    template<std::size_t n>
    Matrix<double,n,n> sampleCovariance(const std::vector<std::vector<double>> &samples){
        int n_dim, m_dim;
        n_dim = samples.size();
        m_dim = samples[0].size();
        MatrixXd data(n_dim,m_dim);
        data = MatVec2EigMat(samples);
        MatrixXd centered = data.rowwise() - data.colwise().mean(); 
        return (centered.adjoint() * centered) / double(data.rows() - 1);
    };  

    template<std::size_t n, std::size_t m>
    double estimateA(const Matrix<double,n,n> &cov){
        SelfAdjointEigenSolver<Matrix<double,n,n>> es;
        es.compute(cov);
        Matrix<double,n,1> eigvals;   
        eigvals = es.eigenvalues().transpose();
        std::sort(eigvals.data(), eigvals.data() + eigvals.size());
        return (eigvals.bottomRows(n-m) * 12).array().sqrt().prod();
    }   

    template<std::size_t n>
    void updateMean(const Matrix<double,n,1>  &sample, std::size_t t, Matrix<double,n,1> &mean){
        mean = mean + (1./(1.+t))*(sample - mean);
    }
    


    template<std::size_t n>
    void updateCovariance(const Matrix<double,n,1> &sample, std::size_t t, Matrix<double,n,1> &mean, Matrix<double,n,n> &cov){
        //double s = (t/((1.+t)*(1.+t)))*(sample-mean).transpose()*(sample-mean); 
        std::size_t t_1 = t+1;
        Matrix<double,n,1> delta_at_nMin1= sample - mean;
        utils::updateMean<n>(sample,t,mean) ;
        Matrix<double,n,1> weighted_delta_at_n = (sample - mean) / t_1;
        Matrix<double,n,n> D_at_n = weighted_delta_at_n.replicate(1,n); 
        Matrix<double,n,n> id;
        id.setIdentity();
        id.diagonal() = delta_at_nMin1;


        //cov = 1./(t+1.) * (cov*(t) + ((t+1.)/(t))*(diff.replicate(1,n))*id);
        cov = cov * ((t_1 - 1.)/ t_1) + (D_at_n*id);
        //cov = (t/(1.+t))*cov; 
        //cov.array() += s;
    }

  /*  template<std::size_t n>
    double squaredDist(const node<n> &n1, const Matrix<double,n,1> &v2){ 
        return squaredDist(n1(), v2);
    };*/
}
// over load add eigen to vector

#endif
