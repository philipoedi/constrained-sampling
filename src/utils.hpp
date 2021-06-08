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
    
    void copyVec2Vec(const std::vector<double>& src, std::vector<double>& dst);

    void copyEig2Vec(const VectorXd& src, std::vector<double>& dst);
    
    std::vector<double> copyEig2Vec(const VectorXd& src);

    void copyMatvec2Matvec(const std::vector<std::vector<double>>& src, std::vector<std::vector<double>>& dst);

    void writeVec2File(const std::vector<std::vector<double>> &data, const std::string &name);
//    void write2file_results(const std::vector<std::vector<double>>& seeds, const std::vector<std::vector>>& samples);
    void writeMetadata2File(const std::string &problem_file, const std::string &name);

    std::string getDateString();
}
// over load add eigen to vector

#endif
