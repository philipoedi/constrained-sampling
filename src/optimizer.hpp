/**
    variety of optimizers (Biased, slack) to use for finding feasible samples
    in a constrained space 
    @file optimizer.hpp
    @author Pihilip Oedi
    @date 2021-06-03
*/

#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <nlopt.hpp>
#include <algorithm>
#include <string>

template<std::size_t n, std::size_t m>
class BaseSampler;

template<std::size_t n, std::size_t m>
class UniformSampler;



#include "sampler.hpp"
#include "utils.hpp"
#include "constraints.hpp"
#include "tangent.hpp"
#include "objectives.hpp"


using namespace nlopt;
using namespace Eigen;


/* todo
    -adjust jacobian for slack opt when using quadratic sonstraint
*/


template<std::size_t n>
class BaseOptimizer
{
    typedef Matrix<double,n,1> Vector;

    public:

        BaseOptimizer();
        BaseOptimizer(const std::vector<double> &lb, const std::vector<double> &ub);
        BaseOptimizer(
            ConstraintCoeffs<n>& eqcons,
            ConstraintCoeffs<n>& ineqcons, 
            const std::vector<double>& lb, 
            const std::vector<double>& ub); 
        BaseOptimizer(
            ConstraintCoeffs<n>& cons,
            const std::vector<double>& lb, 
            const std::vector<double>& ub); 
        void run(const int niter);
        void addConstraints(
            ConstraintCoeffs<n>& eqcons,
            ConstraintCoeffs<n>& ineqcons);
        void addConstraints(std::vector<ConstraintCoeffs<n>> &cons);
        void addConstraints(std::vector<ConstraintCoeffs<n>*> &cons);
        void addConstraints(ConstraintCoeffs<n>& cons);
        void setBounds(const std::vector<double>& lb, const std::vector<double>& ub);
        void results(std::vector<std::vector<double>>& dst);
        std::vector<std::vector<double>> results();
        void samples(std::vector<std::vector<double>>& dst);
        std::vector<std::vector<double>> samples();
        void saveResults(const std::string &name);
        void save_samples(const std::string &name);
        virtual void run(std::vector<std::vector<double>> &seeds){};
        virtual std::vector<double> optimize(std::vector<double> &seed){ return seed;};
        std::vector<int> getNumIterations();

   protected:

        std::vector<int> num_iterations_;
        std::vector<std::vector<double>> results_;
        std::vector<std::vector<double>> samples_;
        opt opt_{"AUGLAG_EQ",n};
        opt local_opt_{"LD_SLSQP",n};
        double minf_;
        UniformSampler<n,0> uni_;
};


template<std::size_t n>
class BiasedOptimizer: public BaseOptimizer<n>
{
    public:

        BiasedOptimizer();
        BiasedOptimizer(
            ConstraintCoeffs<n>& eqcons,
            ConstraintCoeffs<n>& ineqcons, 
            const std::vector<double>& lb, 
            const std::vector<double>& ub): 
            BaseOptimizer<n>(eqcons, ineqcons,lb, ub){}; 
        BiasedOptimizer(
            ConstraintCoeffs<n>& cons,
            const std::vector<double>& lb, 
            const std::vector<double>& ub):
            BaseOptimizer<n>(cons,lb, ub){}; 
        BiasedOptimizer(
            const std::vector<double> &lb,
            const std::vector<double> &ub) : 
            BaseOptimizer<n>(lb, ub){};
        void run(const int niter);
        virtual void run(std::vector<std::vector<double>> &seeds);
        virtual std::vector<double> optimize(std::vector<double> &seed);

    private:

        Bias<n> b_;
};

template<std::size_t n, std::size_t m, std::size_t l>
class SlackOptimizer: public BaseOptimizer<n+m+l+l>
{
    /*
        decision variable vector as follows:

        [x, slack_ineq, slack_eq, t] t ~ introduced variable for l1 norm optimization of slack_eq

    */
    public:

        SlackOptimizer();
        SlackOptimizer(
            ConstraintCoeffs<n>& ineqcons,
            ConstraintCoeffs<n>& eqcons,
            const std::vector<double>& lb,
            const std::vector<double>& ub);
        SlackOptimizer(
            ConstraintCoeffs<n>& cons,
            const std::vector<double>& lb,
            const std::vector<double>& ub);
        void run(const int& niter);
        //void save();
        void addConstraints(ConstraintCoeffs<n>& cons);
        void setBounds(const std::vector<double>& lb, const std::vector<double>& ub);
        double findSlack(const std::vector<double>& x, ConstraintCoeffs<n>& coeffs);
        double findSlack(const std::vector<double>& x, ConstraintCoeffs<n>* coeffs);
        void sample(std::vector<double>& x);

    private:

        SlackData<n,m,l> up_;
        int ineq_count_{0};
        int eq_count_{0};
        std::vector<ConstraintCoeffs<n>*> ineq_cons_;
        std::vector<ConstraintCoeffs<n>*> eq_cons_;
        std::vector<ConstraintCoeffs<n+m+l+l>> ineq_cons_ex_;
        std::vector<ConstraintCoeffs<n+m+l+l>> eq_cons_ex_;
        std::vector<ConstraintCoeffs<n+m+l+l>> eq_cons_ex_up_;
        std::vector<ConstraintCoeffs<n+m+l+l>> eq_cons_ex_low_;
};

// base optimizer member functions
//
//


template<std::size_t n>
BaseOptimizer<n>::BaseOptimizer()
{
    std::cout<<"base"<<std::endl;
    local_opt_.set_xtol_rel(1e-4);
    opt_.set_xtol_rel(1e-4);
    opt_.set_local_optimizer(local_opt_);
    std::cout<<"base2"<<std::endl;
}

template<std::size_t n>
BaseOptimizer<n>::BaseOptimizer(const std::vector<double> &lb, const std::vector<double> &ub) {
    BaseOptimizer();
    setBounds(lb,ub);
};


template<std::size_t n>
BaseOptimizer<n>::BaseOptimizer(
    ConstraintCoeffs<n>& eqcons,
    ConstraintCoeffs<n>& ineqcons,
    const std::vector<double>& lb,
    const std::vector<double>& ub)
    : BaseOptimizer()
{
    this->setBounds(lb, ub);
    this->addConstraints(eqcons, ineqcons);
}

template<std::size_t n>
BaseOptimizer<n>::BaseOptimizer(
    ConstraintCoeffs<n>& cons,
    const std::vector<double>& lb,
    const std::vector<double>& ub)
    : BaseOptimizer()
{
    this->setBounds(lb, ub);
    this->addConstraints(cons);
}

template<std::size_t n>
void BaseOptimizer<n>::setBounds(const std::vector<double>& lb, const std::vector<double>& ub)
{
    uni_.setBounds(lb,ub);
    opt_.set_upper_bounds(ub);
    opt_.set_lower_bounds(lb);
}

template<std::size_t n>
std::vector<int> BaseOptimizer<n>::getNumIterations(){
    return num_iterations_;
}


template<std::size_t n>
void BaseOptimizer<n>::run(const int niter)
{
    std::cout << "base run" << std::endl;   
}

template<std::size_t n>
void BaseOptimizer<n>::addConstraints(
    ConstraintCoeffs<n>& eqcons,
    ConstraintCoeffs<n>& ineqcons)
{
    assert (eqcons.constype == "linear" || eqcons.constype == "quadratic");
    if (eqcons.constype == "linear")
    {
        opt_.add_equality_constraint(linearConstraint<n>, &eqcons, 1e-8);
    }
    else 
    {
        opt_.add_equality_constraint(quadraticConstraint<n>, &eqcons, 1e-8);
    }
    assert (ineqcons.constype == "linear" || ineqcons.constype == "quadratic");
    if (ineqcons.constype == "linear") 
    {
        opt_.add_inequality_constraint(linearConstraint<n>, &ineqcons, 1e-8);
    }
    else
    {
       opt_.add_inequality_constraint(quadraticConstraint<n>, &ineqcons, 1e-8);
    }
}

template<std::size_t n>
void BaseOptimizer<n>::addConstraints(ConstraintCoeffs<n>& cons)
{
    assert (cons.constype == "linear" || cons.constype == "quadratic");
    assert (cons.type == "eq" || cons.type == "ineq");
    if (cons.type == "eq")
    {
        if (cons.constype == "linear")
        {
            opt_.add_equality_constraint(linearConstraint<n>, &cons, 1e-8);
        }
        else 
        {
            opt_.add_equality_constraint(quadraticConstraint<n>, &cons, 1e-8);
        }
    }
    else 
    {
        if (cons.constype == "linear")
        {
            std::cout << "constraint added " << std::endl;
            opt_.add_inequality_constraint(linearConstraint<n>, &cons, 1e-8);
        }
        else 
        {
            opt_.add_inequality_constraint(quadraticConstraint<n>, &cons, 1e-8);
        }
    }
}

template<std::size_t n>
void BaseOptimizer<n>::addConstraints(std::vector<ConstraintCoeffs<n>> &cons){
    for (int i=0; i<cons.size(); i++){
        addConstraints(cons[i]);
    }
}

template<std::size_t n>
void BaseOptimizer<n>::addConstraints(std::vector<ConstraintCoeffs<n>*> &cons){
    for (int i=0; i<cons.size(); i++){
        addConstraints(*(cons[i]));
    }
}

template<std::size_t n>
void BaseOptimizer<n>::results(std::vector<std::vector<double>>& dst)
{
    utils::copyMatvec2Matvec(results_, dst);
}

template<std::size_t n>
std::vector<std::vector<double>> BaseOptimizer<n>::results()
{
    return results_;
}


template<std::size_t n>
void BaseOptimizer<n>::samples(std::vector<std::vector<double>>& dst)
{
    utils::copyMatvec2Matvec(samples_, dst);
}

template<std::size_t n>
std::vector<std::vector<double>> BaseOptimizer<n>::samples()
{
    return samples_;
};  


template<std::size_t n>
void BaseOptimizer<n>::saveResults(const std::string &name)
{
    utils::writeVec2File(results_, name);
};


template<std::size_t n>
void BaseOptimizer<n>::save_samples(const std::string &name)
{
    utils::writeVec2File(samples_, name);
};

// BaseOptimizer member functions finished
// BiasedOptimizer below

template<std::size_t n>
BiasedOptimizer<n>::BiasedOptimizer()
{
    std::cout <<"Biased"  <<std::endl;
}

template<std::size_t n>
void BiasedOptimizer<n>::run(const int niter){
    std::vector<double> x(n);
    this->results_.resize(niter);
    this->samples_.resize(niter);
    for (int i=0; i<niter; ++i)
    {
        b_.x0 = this->uni_.sample();
        utils::copyEig2Vec(b_.x0, x);
        this->opt_.set_min_objective(BiasedObjective<n>, &b_);
        try
        {
            this->opt_.optimize(x, this->minf_);
            this->results_[i].resize(n);
            this->samples_[i].resize(n);
            utils::copyVec2Vec(x, this->results_[i]);
            utils::copyEig2Vec(b_.x0, this->samples_[i]);
        }
        catch(std::exception &e)
        {
            std::cerr << e.what() << std::endl;
        }
        
    }
}

template<std::size_t n>
std::vector<double> BiasedOptimizer<n>::optimize(std::vector<double> &seed){
    b_.x0 = Matrix<double,n,1>(seed.data());
    std::vector<double> x = seed;
    double minf;
    this->opt_.set_xtol_rel(1e-4);
    this->local_opt_.set_xtol_rel(1e-4);
    this->opt_.set_local_optimizer(this->local_opt_);
    this->opt_.set_min_objective(BiasedObjective<n>, &b_);
    try{
        this->opt_.optimize(x,minf);
    }
    catch(std::exception &e){
        std::cerr << e.what() << std::endl;
    }
    this->num_iterations_.push_back(b_.num_iterations);
    b_.num_iterations = 0;
    return x;
}

template<std::size_t n>
void BiasedOptimizer<n>::run(std::vector<std::vector<double>> &seeds){
    for (int i=0; i<seeds.size(); i++){
        (this->results_).push_back(this->optimize(seeds[i]));
    }
}


// slack optimizer below
//
//

template<std::size_t n, std::size_t m, std::size_t l>
SlackOptimizer<n,m,l>::SlackOptimizer(): BaseOptimizer<n+m+l+l>::BaseOptimizer(){};


template<std::size_t n, std::size_t m, std::size_t l>
SlackOptimizer<n,m,l>::SlackOptimizer(
            ConstraintCoeffs<n>& cons,
            const std::vector<double>& lb,
            const std::vector<double>& ub)
            : BaseOptimizer<n+m+l+l>::BaseOptimizer()
{
    this->setBounds(lb, ub);
    this->addConstraints(cons);
};

template<std::size_t n, std::size_t m, std::size_t l>
SlackOptimizer<n,m,l>::SlackOptimizer(
            ConstraintCoeffs<n>& ineqcons,
            ConstraintCoeffs<n>& eqcons,
            const std::vector<double>& lb,
            const std::vector<double>& ub)
            : BaseOptimizer<n+m+l+l>::BaseOptimizer()
{
    this->setBounds(lb, ub);
    this->addConstraints(ineqcons);
    this->addConstraints(eqcons);
};


template<std::size_t n, std::size_t m, std::size_t l>
void SlackOptimizer<n,m,l>::setBounds(const std::vector<double>& lb, const std::vector<double>& ub)
{
    std::vector<double> lb_new;
    std::vector<double> ub_new;
    utils::copyVec2Vec(lb, lb_new);
    utils::copyVec2Vec(ub, ub_new);
    for (int i=0; i<m; i++)
    {
        lb_new.push_back(0);
        ub_new.push_back(HUGE_VAL);
    }
    for (int i=0; i<l+l; i++)
    {
        lb_new.push_back(-HUGE_VAL);
        ub_new.push_back(HUGE_VAL);
    }
    BaseOptimizer<n+m+l+l>::setBounds(lb_new, ub_new);
}

template<std::size_t n, std::size_t m, std::size_t l>
void SlackOptimizer<n,m,l>::addConstraints(ConstraintCoeffs<n>& cons)
{
    assert (cons.type == "ineq" || cons.type == "eq");
    ConstraintCoeffs<n+m+l+l> cons_new;
    
    for (int i=0; i<n; i++){
        cons_new.coeffs(i) = cons.coeffs(i);
    };
    cons_new.cons = cons.cons;
    cons_new.r = cons.r;
    cons_new.type = cons.type;
    cons_new.constype = cons.constype;
    cons_new.sign = cons.sign;

    if (cons.type == "ineq")
    {
        cons_new.coeffs(n+ineq_count_) = -1;
        ineq_count_ += 1;
        ineq_cons_.push_back(&cons);
        ineq_cons_ex_.push_back(cons_new);
        BaseOptimizer<n+m+l+l>::addConstraints(ineq_cons_ex_[ineq_count_-1]);
    }
    else 
    {
        cons_new.coeffs(n+m+eq_count_) = -1;
        ConstraintCoeffs<n+m+l+l> cons_new_eq_up;
        ConstraintCoeffs<n+m+l+l> cons_new_eq_low;
        cons_new_eq_up.type = "ineq";
        cons_new_eq_low.type = "ineq";
        cons_new_eq_up.constype = "linear";
        cons_new_eq_low.constype = "linear";
        cons_new_eq_up.coeffs(n+m+eq_count_) = 1;
        cons_new_eq_up.coeffs(n+m+eq_count_+l) = -1;
        cons_new_eq_low.coeffs(n+m+eq_count_) = -1;
        cons_new_eq_low.coeffs(n+m+eq_count_+l) = -1;
        eq_count_ += 1;
        eq_cons_.push_back(&cons);
        eq_cons_ex_.push_back(cons_new);
        eq_cons_ex_up_.push_back(cons_new_eq_up);
        eq_cons_ex_low_.push_back(cons_new_eq_low);
        BaseOptimizer<n+m+l+l>::addConstraints(eq_cons_ex_[eq_count_-1]);
        BaseOptimizer<n+m+l+l>::addConstraints(eq_cons_ex_up_[eq_count_-1]);
        BaseOptimizer<n+m+l+l>::addConstraints(eq_cons_ex_low_[eq_count_-1]);
    }
}



template<std::size_t n, std::size_t m, std::size_t l>
double SlackOptimizer<n,m,l>::findSlack(const std::vector<double>& x, ConstraintCoeffs<n>& coeffs)
{
    double slack;
    std::vector<double> grad;
    if (coeffs.constype == "linear")
    {
        slack = linearConstraint<n>(x, grad, &coeffs);
    }
    else
    {
        slack = quadraticConstraint<n>(x, grad, &coeffs);
    }
    return slack;
}


template<std::size_t n, std::size_t m, std::size_t l>
double SlackOptimizer<n,m,l>::findSlack(const std::vector<double>& x, ConstraintCoeffs<n>* coeffs)
{
    double slack;
    std::vector<double> grad;
    if (coeffs->constype == "linear")
    {
        slack = linearConstraint<n>(x, grad, coeffs);
    }
    else
    {
        slack = quadraticConstraint<n>(x, grad, coeffs);
    }
    return slack;
}


template<std::size_t n, std::size_t m, std::size_t l>
void SlackOptimizer<n,m,l>::sample(std::vector<double> &x)
{
    // loop over ineq constraint coeffs find slack
    ConstraintCoeffs<n>*  coeffs_temp;
    double slack{0}, slack_temp{0};
    std::vector<double> x0;
    utils::copyEig2Vec(this->uni_.sample(),x);
    utils::copyVec2Vec(x, x0);
    x0.resize(n);
    //utils::copyEig2Vec(this->uni_.sample().head(n),x0);
    // initialize slack for inequalities
    for (int i=0; i<m; i++)
    {
        coeffs_temp = ineq_cons_[i];
        slack_temp = this->findSlack(x0 ,coeffs_temp);
        slack_temp = (slack_temp < 0) ? 0 : slack_temp;
        slack = (slack_temp > slack) ? slack_temp : slack;
        x[n+i] = slack;
            
    };
    slack = 0;
    for (int i=0; i<l; i++)
    {
        coeffs_temp = eq_cons_[i];
        slack = this->findSlack(x0 ,coeffs_temp);
        x[n+m+i] = slack;
        // choose t_i to be abs(slack_i)
        x[n+m+l+i] = (slack < 0) ? slack*-1 : slack;
    };
}


template<std::size_t n, std::size_t m, std::size_t l>
void SlackOptimizer<n,m,l>::run(const int& niter)
{
    initSlackData<n,m,l>(up_);
    std::vector<double> x(n+m+l+l);
    this->results_.resize(niter);
    this->samples_.resize(niter);
    for (int i=0; i<niter; ++i)
    {
        this->sample(x);
        this->samples_[i].resize(n+m+l+l);
        utils::copyVec2Vec(x, this->samples_[i]);
        this->opt_.set_min_objective(slackObjective<n,m,l>, &up_);
        try
        {
            this->opt_.optimize(x, this->minf_);
            this->results_[i].resize(n);
            utils::copyVec2Vec(x, this->results_[i]);
        }
        catch(std::exception &e)
        {
            std::cerr << e.what() << std::endl;
        }
    }
}

#endif



