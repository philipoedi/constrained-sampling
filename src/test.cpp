#include <Eigen/Dense>
#include <vector>
#include <iostream>

using namespace Eigen;

int main()
{
    std::vector<double> a{1,2};
    Vector2d b;
    for (std::size_t i = 0; i < a.size(); ++i){
        b(i) = a[i];
    };
    std::cout << b << std::endl;

}
