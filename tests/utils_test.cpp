#include <iostream>
#include "utils.hpp"
#include <cassert>
#include <vector>
#include <algorithm>
#include <functional>



int main(){
    using namespace std::placeholders;    

    // test factorial
    const int n = 3;
    int res;
    res = utils::factorial(n);
    std::cout << res << std::endl;
    assert (res == 6);
    

    std::vector<double> vec1{1,2,3};
    std::vector<double> vec2{2,3,4};
    utils::copyVec2Vec(vec1, vec2);
    for (int i=0; i<vec2.size(); i++)
    {
        std::cout << vec2[i] << std::endl;
        assert (vec1[i] == vec2[i]);
    }

    std::vector<double> vec3{2,3,4};
    Vector3d eig1(vec3.data());
    Vector3d eig2(vec1.data());
    utils::copyEig2Vec(eig1, vec1);
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
    utils::copyMatvec2Matvec(mat1, mat2);
     for (int i=0; i<mat1.size(); i++)
    {
        assert (mat1[i][0] = mat2[i][0]);
        assert (mat1[i][1] = mat2[i][1]);
    }
     
    // test writeVec2File
    std::cout << "test writing to file" << std::endl;
    std::vector<std::vector<double>> data(2,std::vector<double>(4));
    std::string name{"example_file"};
    //utils::writeVec2File<double>(data, name);
    
    name = "../src/linear_cons.cpp";
    utils::writeMetadata2File(name,"example_file");
    
    std::string today = utils::getDateString();
    std::cout << today << std::endl;
    
    std::string todaytime = utils::getDateTimeString();
    std::cout << todaytime << std::endl;
    std::cout << "testing copy_eig2vec" << std::endl;

    std::vector<double> vec = utils::copyEig2Vec(eig1);
    std::vector<double>::iterator it;
    for (it = vec.begin(); it != vec.end(); ++it)
    {
        std::cout << *it << std::endl;
    }

    std::cout << "$distance: "<< utils::squaredDist<n>(eig1, eig2) << std::endl;
    assert(utils::squaredDist<n>(eig1,eig2) == 3);
    std::vector<Vector3d> eigs;
    eigs.push_back(eig1);
    eigs.push_back(eig1);
    eigs.push_back(eig1);
    std::vector<double> out;
    std::transform(eigs.begin(),eigs.end(),std::back_inserter(out),std::bind(utils::squaredDist<n>, eig2, _1));
    std::cout << out[0] << out[2] << std::endl;

    double ress[3];
    std::cout << "testing copyEig2Arr" << std::endl;
    utils::copyEig2Arr(eig2,ress);
    std::cout << ress[0] << " " << ress[1] << " " << ress[2] << std::endl;
    std::cout << eig2 <<std::endl;
    
    std::vector<double> xx(3,1);
    std::vector<double> uu = utils::slice(xx,0,1);    
    std::cout <<"size of uu: "<<  uu.size() << std::endl;

    return 0;
}

