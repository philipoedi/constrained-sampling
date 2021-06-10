#include <vector>
#include <Eigen/Dense>
#include "RRT.hpp"


int main(){
    using namespace Eigen;
    const std::size_t n{2};
    std::vector<double> vec1{1,2};
    Vector2d eig1(vec1.data()); 
    node<n> n1;
    n1.location = eig1;
    std::cout << n1() << std::endl; 

    return 0;
};
