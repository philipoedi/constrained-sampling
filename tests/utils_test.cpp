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
    return 0;
}

