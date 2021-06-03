/**
    variety of optimizers (Biased, slack) to use for finding feasible samples
    in a constrained space 
    @file optimizer.hpp
    @author Pihilip Oedi
    @date 2021-06-03
*/

#ifndef OPTIMIZER_H
#define OPTIMIZER_H
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <nlopt.hpp>
#include <list>
#include <string>
#include "sampler.hpp"
#include "utils.hpp"

using namespace nlopt;
using namespace Eigen;


/* todo
    -adjust jacobian for slack opt when using quadratic sonstraint
*/


/**
    Specifies the vector of Biases in Biased optimization.
    @tparam n Dimension N of vector
*/
template<std::size_t n>
struct Bias {
    /// Eigen Vector representing the bias
    Matrix<double, n, 1> x0;
};

/**
    Auxilliary data used to specifie constraints, such as coefficients,
    constraint_type("ineq" or "eq") and type "linear" or "quadratic"
    @tparam n Dimension N of vector
*/
template<std::size_t n>
struct ConstraintCoeffs {
    /// 
    Matrix<double, n, 1> coeffs = Matrix<double, n, 1>::Zero();
    Matrix<double, n, 1> q = Matrix<double, n, 1>::Zero();
    Matrix<double, n, n> P = Matrix<double, n, n>::Zero();
    double cons;
    double sign;
    double r{0};
    std::string type;
    std::string constype;
};

template<std::size_t n, std::size_t m, std::size_t l>
struct SlackData {
    bool state{false};
    Matrix<double, n+m+l+l, 1> a = Matrix<double, n+m+l+l,1>::Zero();
    Matrix<double, n+m+l+l, 1> g = Matrix<double, n+m+l+l,1>::Zero();
};

/**
    Initializes slack variables to use for slack optimization
    @param sl SlackData struct containing vectors of decision variables
*/
template<std::size_t n, std::size_t m, std::size_t l>
void initSlackData(SlackData<n,m,l>& sl)
{
    for (int i=0; i<n+m+l+l; i++)
    {
        if (i >= n && i <n+m)
        {
            sl.a(i) = 1;
            sl.g(i) = 1;
        }
        else if (i >= n+m+l)
        {
            sl.a(i) = 1;
            sl.g(i) = 1;
        };
    };
};

/**
    Evaluates objective function in Biased optimization
    @param x Current location at which to evaluate objective
    @param grad Vector of current gradient
    @param data Auxilliary data to be used in calculations, of type Bias
    @return Objective value
*/
template<std::size_t n> 
double BiasedObjective(const std::vector<double>& x, std::vector<double>& grad, void*data )
{
    typedef Matrix<double, n, 1> vec;
    Bias<n> *b = (Bias<n>*) data;
    vec x_vec(x.data());
    vec x_x0 = x_vec - b->x0;
    if(!grad.empty()){       
        utils::copyEig2Vec(2*x_x0, grad);
/*    for (std::size_t i = 0; i<x_x0.size() ;++i){
        grad[i] = 2*x_x0[i];
    }*/
    }   
    return x_x0.transpose()*x_x0;    
}

/**
    Evaluates objective function in slack optimization
    @param x Current location at which to evaluate objective
    @param grad Vector of current gradient
    @param data Auxilliary data to be used in calculations
    @return Objective value
*/
template<std::size_t n, std::size_t m, std::size_t l> 
double slackObjective(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    typedef Matrix<double, n+m+l+l,1> vec;
    SlackData<n,m,l>* u = (SlackData<n,m,l>*) data;
    vec x_vec(x.data());
    if(!grad.empty()){
        utils::copyEig2Vec(u->a, grad);
    }
    return u->a.transpose()*x_vec;
}

/**
    Evaluates a linear constraint, function can be used with nlopt
    @param x Current location at which to evaluate constraint 
    @param grad Vector of current gradient
    @param data Auxilliary data to be used in calculations
    @return Constraint value
*/
template<std::size_t n> 
double linearConstraint(const std::vector<double>& x, std::vector<double> &grad, void*data)
{
    typedef Matrix<double, n, 1> vec;
    ConstraintCoeffs<n> *c = (ConstraintCoeffs<n>*) data;
    vec x_vec(x.data());
    if (!grad.empty()){
        utils::copyEig2Vec(c->coeffs, grad);
    }
    return x_vec.transpose() * c->coeffs - c->cons;
}

/**
    Evaluates a linear constraint, function can be used with nlopt
    @param x Current location at which to evaluate constraint 
    @param grad Vector of current gradient
    @param data Auxilliary data to be used in calculations
    @return Constraint value
*/
template<std::size_t n> 
double linearConstraint(const Matrix<double,n,1> &x, Matrix<double,n,1> &grad, void*data)
{
    typedef Matrix<double, n, 1> vec;
    ConstraintCoeffs<n> *c = (ConstraintCoeffs<n>*) data;
    if (!grad.empty()){
        grad = c->coeffs;
        //utils::copyEig2Vec(c->coeffs, grad);
    }
    return x.transpose() * c->coeffs - c->cons;
}

/**
    Evaluates a linear constraint, function can be used with nlopt
    @param x Current location at which to evaluate constraint 
    @param grad Vector of current gradient
    @param data Auxilliary data to be used in calculations
    @return Constraint value
*/
template<std::size_t n> 
double linearConstraint(const Matrix<double,n,1> &x, void*data)
{
    ConstraintCoeffs<n> *c = (ConstraintCoeffs<n>*) data;
    return x.transpose() * c->coeffs - c->cons;
}

/**
    Evaluates a quadratic constraint, function can be used with nlopt
    @param x Current location at which to evaluate constraint 
    @param grad Vector of current gradient
    @param data Auxilliary data to be used in calculations, of type ConstraintCoeffs that
        specifies coefficients of constraint
    @return Constraint value
*/
template<std::size_t n> 
double quadraticConstraint(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    typedef Matrix<double, n, 1> vec;
    double res = 0;
    ConstraintCoeffs<n> *c = (ConstraintCoeffs<n>*) data;
    vec x_vec(x.data());
    // 0.5 * x.T@P@x + q.T@x+ r
    if (!grad.empty()){
        utils::copyEig2Vec(c->P.transpose()*x_vec + c->q, grad);
    }
    res += 0.5 * x_vec.transpose() * c->P * x_vec; 
    res += x_vec.transpose()*c->q;  
    return res - c->r;
}

template<std::size_t n>
class BaseOptimizer
{
    typedef Matrix<double,n,1> Vector;

    public:

        BaseOptimizer();
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
        void addConstraints(ConstraintCoeffs<n>& cons);
        void setBounds(const std::vector<double>& lb, const std::vector<double>& ub);
        void results(std::vector<std::vector<double>>& dst);
        std::vector<std::vector<double>> results();
        void samples(std::vector<std::vector<double>>& dst);
        std::vector<std::vector<double>> samples();
        void saveResults(const std::string &name);
        void save_samples(const std::string &name);

   protected:

        std::vector<std::vector<double>> results_;
        std::vector<std::vector<double>> samples_;
        opt opt_{"AUGLAG_EQ",n};
        opt local_opt_{"LD_SLSQP",n};
        double minf_;
        UniformSampler<n> uni_;
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
      void run(const int niter);
    
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
}

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



