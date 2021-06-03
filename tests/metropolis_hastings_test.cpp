#include <functional>
#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <random>
#include "sampler.hpp"
#include "optimizer.hpp"
#include <list>

using namespace Eigen;


const std::size_t n{2};
const int n_iter{5};

typedef Matrix<double,n,1> Vector;


double q(const Vector2d &x, const Vector2d &x2)
{
    return 1;
}


template<std::size_t n>
struct p
{   
    // evaluate constraints, if constraints not satisfied return prob density of 0
    // else returns 1
    double slack;
    std::vector<ConstraintCoeffs<n>> cons;
    double operator()(const Matrix<double, n, 1> &x)
    {
        std::vector<double> x_vec;
        utils::copyEig2Vec(x, x_vec);
        typename std::vector<ConstraintCoeffs<n>>::iterator it;
        for (it = cons.begin(); it != cons.end(); ++it)
        {
            slack = linearConstraint<n>(x_vec, *it);
            if (slack > 0)
                return 0;
        }
        return 1;
    };
};


template<std::size_t n>
struct proposal_specs
{
    UniformSampler<n> uni;
    Matrix<double,n,1> lb, ub, lb_new, ub_new, widths;
    Matrix<double,n,1> operator()(const Matrix<double,n,1> &x)
    {
       lb_new = x - widths;
       ub_new = x + widths;
       lb_new = lb_new.cwiseMax(lb);
       ub_new = ub_new.cwiseMin(ub);
       uni.setBounds(lb_new, ub_new); 
       return uni.sample();
    };
};

int main()
{
    std::vector<double> lb{0,0};
    std::vector<double> ub{2,2};
   
//    uniform_neighborhood_sampler proposal_sampler(lb, ub);
  //  Vector (*Q)(const Vector&) = &proposal_sampler::sample;
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0,3);
    std::cout << "generating 20 random numbers between 0 and 3" << std::endl;
    for (int i=0; i< 20; i++)
        std::cout << dist(gen) << std::endl;
    Vector2d lb_eig(lb.data()); 
    Vector2d ub_eig(ub.data()); 
    Vector2d x;
    x << 1.,  1.;
    Vector2d x2{0,0};
    Vector2d w{0.5,0.5};

    std::cout << "default lower bounds of ps" << std::endl;
    proposal_specs<n> ps;
    ps.lb = lb_eig;
    ps.ub = ub_eig;
    ps.widths = w;
    std::cout << ps.lb << std::endl;
    std::cout << "trying functional of proposal_specs"<< std::endl;
    std::cout << ps(x)<< std::endl;

    std::cout << "x2 = (0,0)" << std::endl;
    for (int i=0; i< 20; i++)
        std::cout << ps(x2) << std::endl;

    std::cout << "initializing metropolis hastings" << std::endl; 
    MetropolisHastings<n> mh(lb, ub); 
    std::function<Matrix<double,n,1>(const Matrix<double,n,1>&)> f = ps;
    std::cout << f(x) << std::endl;
    mh.set_Q(f);
    std::cout << "trying Q usign std::function" << std::endl;
    std::cout << mh.sampleQ(x) << std::endl;
    std::cout << "evaluate constraints and prob density of true dist" << std::endl;
    ConstraintCoeffs<n> cons;
    p<n> p_dense;
    mh.setP(p_dense);
    mh.setQ(q);
    std::cout << "trying run" << std::endl;
    mh.run(n_iter); 
    return 0;
}
