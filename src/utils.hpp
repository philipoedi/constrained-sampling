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
   

  /*  template<std::size_t n>
    double squaredDist(const node<n> &n1, const Matrix<double,n,1> &v2){ 
        return squaredDist(n1(), v2);
    };*/
}
// over load add eigen to vector

#endif
