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
        //assert (src.size() == dst.size());
        dst.resize(src.size());
        for (int i=0; i<src.size(); i++)
        {
            dst[i] = src[i];
        }
    };

    void copy_eig2vec(const VectorXd& src, std::vector<double>& dst)
    {
        assert (src.size() == dst.size());
        for (int i=0; i<src.size(); i++)
        {
            dst[i] = src(i);
        }
    };
    
    void copy_matvec2matvec(const std::vector<std::vector<double>>& src, std::vector<std::vector<double>>& dst)
    {
        std::size_t d{src.size()};
        dst.resize(d);
        for (int i=0; i<d; i++)
        {
            dst[i].resize(src[i].size());
            copy_vec2vec(src[i], dst[i]);
        }
    };
}
