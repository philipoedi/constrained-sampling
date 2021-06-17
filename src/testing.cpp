#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;
int main(){
    double foo[6]= {1,2,3,4,5,6}; 
    Matrix<double,2,3>k(foo); 
    cout << Matrix<double,2,3>(foo) << endl;
    foo = (k*2).data();
    std::cout << foo[0] << std::endl; 
    return 0;
};
