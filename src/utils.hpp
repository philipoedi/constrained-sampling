#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <vector>


using namespace Eigen;
/*
string combine_sample_result();

void write_to_file(ofstream file, string results);


void run_sampling();
*/
namespace utils
{
    int factorial(int n);
    
    void copy_vec2vec(const std::vector<double>& src, std::vector<double>& dst);

    void copy_eig2vec(const VectorXd& src, std::vector<double>& dst);

//    void write2file_results(const std::vector<std::vector<double>>& seeds, const std::vector<std::vector>>& samples);
}
// over load add eigen to vector

#endif
