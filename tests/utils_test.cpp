#include <iostream>
#include "utils.hpp"
#include <cassert>
#include <vector>


int main(){
    

    // test factorial
    int n = 3;
    int res;
    res = utils::factorial(n);
    std::cout << res << std::endl;
    assert (res == 6);
    

    std::vector<double> vec1{1,2,3};
    std::vector<double> vec2{2,3,4};
    utils::copy_vec2vec(vec1, vec2);
    for (int i=0; i<vec2.size(); i++)
    {
        std::cout << vec2[i] << std::endl;
        assert (vec1[i] == vec2[i]);
    }

    std::vector<double> vec3{2,3,4};
    Vector3d eig1(vec3.data());
    utils::copy_eig2vec(eig1, vec1);
    for (int i=0; i<vec1.size(); i++)
     {
        std::cout << vec1[i] << std::endl;
        assert (vec1[i] == eig1(i));
     }
    std::vector<std::vector<double>> mat1;
    std::vector<std::vector<double>> mat2;
    mat1.resize(2);
    for (int i=0; i<mat1.size(); i++)
    {
        mat1[i].resize(2);
        mat1[i][0] = 1;
        mat1[i][1] = 2;
    }
    utils::copy_matvec2matvec(mat1, mat2);
     for (int i=0; i<mat1.size(); i++)
    {
        assert (mat1[i][0] = mat2[i][0]);
        assert (mat1[i][1] = mat2[i][1]);
    }
     
    return 0;
}

