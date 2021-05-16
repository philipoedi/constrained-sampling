#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <vector>
#include <string>


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
    
    void copy_matvec2matvec(const std::vector<std::vector<double>>& src, std::vector<std::vector<double>>& dst);

    void write_vec2file(const std::vector<std::vector<double>> &data, const std::string &name);
//    void write2file_results(const std::vector<std::vector<double>>& seeds, const std::vector<std::vector>>& samples);
    void write_metadata2file(const std::string &problem_file, const std::string &name);

    std::string get_date_string();
}
// over load add eigen to vector

#endif
