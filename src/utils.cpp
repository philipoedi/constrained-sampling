#include "utils.hpp"
#include <cassert>
#include <Eigen/Dense>

using namespace Eigen;


namespace utils 
{
    int factorial(int n)
    {
        int result = 1.;
        if (n != 0){
            for (int i= 1; i <= n; ++i){
                result *= i;
            }
        }
        return result;
    };

    void copy_vec2vec(const std::vector<double>& src, std::vector<double>& dst)
    {
        assert (src.size() == dst.size());
        for (int i=0; i<src.size(); i++)
        {
            dst[i] = src[i];
        }
    }

    void copy_eig2vec(const VectorXd& src, std::vector<double>& dst)
    {
        assert (src.size() == dst.size());
        for (int i=0; i<src.size(); i++)
        {
            dst[i] = src(i);
        }
    }

}
