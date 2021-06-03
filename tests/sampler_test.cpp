#include <iostream>
#include "sampler.hpp"
#include <vector>
#include <cassert>
#include <algorithm>


void check_bounds(std::pair<Vector2d,Vector2d> b)
{
    assert (b.first(0) == 0);
    assert (b.second(1) == 2);
}

double func(const Vector2d &a)
{
    return a.sum();
}

double func2(const Vector2d &a, const Vector2d &b)
{
    return a.sum()  + b.sum();
}
  

int f()
{
    int i = 1;
    return i;
}

double P(const Matrix<double,2,1> &x)
{
    return 1;
};


int main()
{
    const std::size_t n = 2;

    // check default init
    UniformSampler<n> uni;
    
    // check setBounds using std::vector
    std::vector<double> ub{1,2};
    std::vector<double> lb{0,0};
    uni.setBounds(lb, ub);
    std::pair<Vector2d,Vector2d> bounds;
    bounds = uni.getBounds();
   // assert (bounds.first(0) == 0);
   // assert (bounds.second(1) == 2);
    check_bounds(bounds);


    // check setBounds using eigen vec
    Vector2d lb_eig(lb.data());
    Vector2d ub_eig(ub.data());
    std::cout << lb_eig << std::endl;
    uni.setBounds(lb_eig, ub_eig);
    bounds = uni.getBounds();
    check_bounds(bounds);
    assert (bounds.first(0) == 0);
    assert (bounds.second(1) == 2);

    // check sample
    Vector2d sample;
    for (int i=0; i<10; i++)
    {
        sample = uni.sample();
        std::cout << "sample: " << sample << std::endl;
    }
    assert (sample(0) <= ub_eig(0));
    assert (sample(1) <= ub_eig(1));
    assert (sample(0) >= lb_eig(0));
    assert (sample(1) >= lb_eig(1));

    MetropolisHastings<n> mh;
    mh.setBounds(lb, ub);
    bounds = mh.getBounds();
    check_bounds(bounds);

    MetropolisHastings<n> mh2(lb, ub);
    bounds = mh.getBounds();
    check_bounds(bounds);
    mh.setP(func);
    mh.setQ(func2);
    MetropolisHastings<n> mh3(lb_eig, ub_eig);
    bounds = mh.getBounds();
    check_bounds(bounds);



    std::cout << mh.getP(ub_eig) << std::endl;
    std::cout << mh.getQ(ub_eig,lb_eig) << std::endl;
    std::cout << "Default A: "<< mh.aDefault(0.22, ub_eig, lb_eig) << std::endl;
    
    int j{0};
    j = f();
    std::cout << j << std::endl;
    /*uniform_neighborhood_sampler<n> uns(lb_eig, ub_eig);
    bounds = uns.getBounds();
    check_bounds(bounds);
    uns.set_P(P);
    std::cout << "checking get_P:" <<  uns.get_P(lb_eig) << std::endl;
    */
    return 0;
}

