#include <vector>
#include <Eigen/Dense>
#include "RRT.hpp"
#include "optimizer.hpp"

int main(){
    using namespace Eigen;
    const std::size_t n{3};
    const std::size_t m{1};
    std::vector<double> vec1{1,2,3};
    Vector3d eig1(vec1.data()); 
    node<n> n1;
    n1.location = eig1;
    std::cout << n1() << std::endl; 
    
    // bounds
    std::vector<double> lb{-2,-2,-2};
    std::vector<double> ub{2,2,2};
    
    // step size
    double alpha{0.5};

    // create RRT obj
    RRT<n,m> rrt;
    rrt.setBounds(lb, ub);
    rrt.setStepSize(alpha);

    // add constraints
    ConstraintCoeffs<n> c;
    c.type = "eq";
    c.constype = "linear";
    std::vector<double> coeff_vec{-1,-1,-1};
    Matrix<double,n,1> coeff(coeff_vec.data());
    c.coeffs = coeff;
    c.cons = -1;         
    std::vector<ConstraintCoeffs<n>> cons;
    cons.push_back(c);
    rrt.addConstraints(cons);
    // run
    const int n_iter{100};
    rrt.run(n_iter);
    // save files
    rrt.saveResults("results");
    rrt.saveSamples("samples"); 

    std::vector<double> seed{0,0,1};
    //rrt.runOnTangent(10, seed, lb, ub); 
    
    RRT<n,m> rrt22;
    rrt22.setStepSize(1e-2);
    BaseSampler<n> * rrt2 = &rrt22; 
    rrt2->setBounds(lb,ub);
    ConstraintCoeffs<n> c2;
    c2.type = "eq";
    c2.constype = "quadratic";
    Matrix<double,n,n> ind = Matrix<double,n,n>::Identity(); 
    c2.P = ind*2;
    c2.r = 1;         
    std::vector<ConstraintCoeffs<n>> cons2;
    cons2.push_back(c2);
    rrt2->addConstraints(cons2);
    BaseOptimizer<n> * b = new BiasedOptimizer<n>(lb, ub);
    rrt2->setOptimizer(b);
    std::cout << rrt2->hasOptimizer() << std::endl;
    rrt22.runOnTangent(10,seed,lb,ub);
    rrt22.setUseTangent(true);
    rrt22.run(10,seed,lb,ub);
   std::cout <<"finnisef2" << std::endl;   // run


    return 0;
};
