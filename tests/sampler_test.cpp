#include <iostream>
#include "sampler.hpp"
#include <vector>
#include <cassert>
#include <algorithm>


int main()
{
    const std::size_t n = 2;

    // check default init
    uniform_sampler<n> uni;
    
    // check set_bounds using std::vector
    std::vector<double> ub{1,2};
    std::vector<double> lb{0,0};
    uni.set_bounds(lb, ub);
    std::pair<Vector2d,Vector2d> bounds;
    bounds = uni.get_bounds();
    assert (bounds.first(0) == 0);
    assert (bounds.second(1) == 2);

    // check set_bounds using eigen vec
    Vector2d lb_eig(lb.data());
    Vector2d ub_eig(ub.data());
    uni.set_bounds(lb_eig, ub_eig);
    bounds = uni.get_bounds();
    assert (bounds.first(0) == 0);
    assert (bounds.second(1) == 2);

    // check sample
    Vector2d sample;
    sample = uni.sample();
    assert (sample(0) <= ub_eig(0));
    assert (sample(1) <= ub_eig(1));
    assert (sample(0) >= lb_eig(0));
    assert (sample(1) >= lb_eig(1));
}
