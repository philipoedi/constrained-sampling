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


template<std::size_t n>
struct bias {
    Matrix<double, n, 1> x0;
};

template<std::size_t n>
struct constraint_coeffs {
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
struct slack_data {
    bool state{false};
    Matrix<double, n+m+l+l, 1> a = Matrix<double, n+m+l+l,1>::Zero();
    Matrix<double, n+m+l+l, 1> g = Matrix<double, n+m+l+l,1>::Zero();
};

template<std::size_t n, std::size_t m, std::size_t l>
void init_slack_data(slack_data<n,m,l>& sl)
// sets the a vector according to the objective function
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

template<std::size_t n> double biased_objective(const std::vector<double>& x, std::vector<double>& grad, void*data )
{
    typedef Matrix<double, n, 1> vec;
    bias<n> *b = (bias<n>*) data;
    vec x_vec(x.data());
    vec x_x0 = x_vec - b->x0;
    if(!grad.empty()){       
        utils::copy_eig2vec(2*x_x0, grad);
/*    for (std::size_t i = 0; i<x_x0.size() ;++i){
        grad[i] = 2*x_x0[i];
    }*/
    }   
    return x_x0.transpose()*x_x0;    
}

template<std::size_t n, std::size_t m, std::size_t l> 
double slack_objective(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    typedef Matrix<double, n+m+l+l,1> vec;
    slack_data<n,m,l>* u = (slack_data<n,m,l>*) data;
    vec x_vec(x.data());
    if(!grad.empty()){
        utils::copy_eig2vec(u->a, grad);
    }
    return u->a.transpose()*x_vec;
}

template<std::size_t n> double linear_constraint(const std::vector<double>& x, std::vector<double> &grad, void*data)
{
    typedef Matrix<double, n, 1> vec;
    constraint_coeffs<n> *c = (constraint_coeffs<n>*) data;
    vec x_vec(x.data());
    if (!grad.empty()){
        utils::copy_eig2vec(c->coeffs, grad);
    }
    return x_vec.transpose() * c->coeffs - c->cons;
}


template<std::size_t n> double quadratic_constraint(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    typedef Matrix<double, n, 1> vec;
    double res = 0;
    constraint_coeffs<n> *c = (constraint_coeffs<n>*) data;
    vec x_vec(x.data());
    // 0.5 * x.T@P@x + q.T@x+ r
    if (!grad.empty()){
        utils::copy_eig2vec(c->P.transpose()*x_vec + c->q, grad);
    }
    res += 0.5 * x_vec.transpose() * c->P * x_vec; 
    res += x_vec.transpose()*c->q;  
    return res - c->r;
}

template<std::size_t n>
class base_optimizer
{
    typedef Matrix<double,n,1> Vector;

    public:

        base_optimizer();
        base_optimizer(
            constraint_coeffs<n>& eqcons,
            constraint_coeffs<n>& ineqcons, 
            const std::vector<double>& lb, 
            const std::vector<double>& ub); 
        base_optimizer(
            constraint_coeffs<n>& cons,
            const std::vector<double>& lb, 
            const std::vector<double>& ub); 
        void run(const int niter);
//        virtual void save();
        void add_constraints(
            constraint_coeffs<n>& eqcons,
            constraint_coeffs<n>& ineqcons);
        void add_constraints(constraint_coeffs<n>& cons);
        void set_bounds(const std::vector<double>& lb, const std::vector<double>& ub);
        void results(std::vector<std::vector<double>>& dst);
        std::vector<std::vector<double>> results();
        void samples(std::vector<std::vector<double>>& dst);
        std::vector<std::vector<double>> samples();
        void save_results(const std::string &name);
        void save_samples(const std::string &name);

   protected:

        std::vector<std::vector<double>> results_;
        std::vector<std::vector<double>> samples_;
        opt opt_{"AUGLAG_EQ",n};
        opt local_opt_{"LD_SLSQP",n};
        double minf_;
        uniform_sampler<n> uni_;
};


template<std::size_t n>
class biased_optimizer: public base_optimizer<n>
{
    public:

        biased_optimizer();
        biased_optimizer(
            constraint_coeffs<n>& eqcons,
            constraint_coeffs<n>& ineqcons, 
            const std::vector<double>& lb, 
            const std::vector<double>& ub): 
            base_optimizer<n>(eqcons, ineqcons,lb, ub){}; 
        biased_optimizer(
            constraint_coeffs<n>& cons,
            const std::vector<double>& lb, 
            const std::vector<double>& ub):
            base_optimizer<n>(cons,lb, ub){}; 
      void run(const int niter);
    
    private:

        bias<n> b_;
};

template<std::size_t n, std::size_t m, std::size_t l>
class slack_optimizer: public base_optimizer<n+m+l+l>
{
    /*
        decision variable vector as follows:

        [x, slack_ineq, slack_eq, t] t ~ introduced variable for l1 norm optimization of slack_eq

    */
    public:

        slack_optimizer();
        slack_optimizer(
            constraint_coeffs<n>& ineqcons,
            constraint_coeffs<n>& eqcons,
            const std::vector<double>& lb,
            const std::vector<double>& ub);
        slack_optimizer(
            constraint_coeffs<n>& cons,
            const std::vector<double>& lb,
            const std::vector<double>& ub);
        void run(const int& niter);
        //void save();
        void add_constraints(constraint_coeffs<n>& cons);
        void set_bounds(const std::vector<double>& lb, const std::vector<double>& ub);
        double find_slack(const std::vector<double>& x, constraint_coeffs<n>& coeffs);
        double find_slack(const std::vector<double>& x, constraint_coeffs<n>* coeffs);
        void sample(std::vector<double>& x);

    private:

        slack_data<n,m,l> up_;
        int ineq_count_{0};
        int eq_count_{0};
        std::vector<constraint_coeffs<n>*> ineq_cons_;
        std::vector<constraint_coeffs<n>*> eq_cons_;
        std::vector<constraint_coeffs<n+m+l+l>> ineq_cons_ex_;
        std::vector<constraint_coeffs<n+m+l+l>> eq_cons_ex_;
        std::vector<constraint_coeffs<n+m+l+l>> eq_cons_ex_up_;
        std::vector<constraint_coeffs<n+m+l+l>> eq_cons_ex_low_;
};

// base optimizer member functions
//
//


template<std::size_t n>
base_optimizer<n>::base_optimizer()
{
    std::cout<<"base"<<std::endl;
    local_opt_.set_xtol_rel(1e-4);
    opt_.set_xtol_rel(1e-4);
    opt_.set_local_optimizer(local_opt_);
}

template<std::size_t n>
base_optimizer<n>::base_optimizer(
    constraint_coeffs<n>& eqcons,
    constraint_coeffs<n>& ineqcons,
    const std::vector<double>& lb,
    const std::vector<double>& ub)
    : base_optimizer()
{
    this->set_bounds(lb, ub);
    this->add_constraints(eqcons, ineqcons);
}

template<std::size_t n>
base_optimizer<n>::base_optimizer(
    constraint_coeffs<n>& cons,
    const std::vector<double>& lb,
    const std::vector<double>& ub)
    : base_optimizer()
{
    this->set_bounds(lb, ub);
    this->add_constraints(cons);
}

template<std::size_t n>
void base_optimizer<n>::set_bounds(const std::vector<double>& lb, const std::vector<double>& ub)
{
    uni_.set_bounds(lb,ub);
    opt_.set_upper_bounds(ub);
    opt_.set_lower_bounds(lb);
}

template<std::size_t n>
void base_optimizer<n>::run(const int niter)
{
    std::cout << "base run" << std::endl;   
}


template<std::size_t n>
void base_optimizer<n>::add_constraints(
    constraint_coeffs<n>& eqcons,
    constraint_coeffs<n>& ineqcons)
{
    assert (eqcons.constype == "linear" || eqcons.constype == "quadratic");
    if (eqcons.constype == "linear")
    {
        opt_.add_equality_constraint(linear_constraint<n>, &eqcons, 1e-8);
    }
    else 
    {
        opt_.add_equality_constraint(quadratic_constraint<n>, &eqcons, 1e-8);
    }
    assert (ineqcons.constype == "linear" || ineqcons.constype == "quadratic");
    if (ineqcons.constype == "linear") 
    {
        opt_.add_inequality_constraint(linear_constraint<n>, &ineqcons, 1e-8);
    }
    else
    {
       opt_.add_inequality_constraint(quadratic_constraint<n>, &ineqcons, 1e-8);
    }
}

template<std::size_t n>
void base_optimizer<n>::add_constraints(constraint_coeffs<n>& cons)
{
    assert (cons.constype == "linear" || cons.constype == "quadratic");
    assert (cons.type == "eq" || cons.type == "ineq");
    if (cons.type == "eq")
    {
        if (cons.constype == "linear")
        {
            opt_.add_equality_constraint(linear_constraint<n>, &cons, 1e-8);
        }
        else 
        {
            opt_.add_equality_constraint(quadratic_constraint<n>, &cons, 1e-8);
        }
    }
    else 
    {
        if (cons.constype == "linear")
        {
            std::cout << "constraint added " << std::endl;
            opt_.add_inequality_constraint(linear_constraint<n>, &cons, 1e-8);
        }
        else 
        {
            opt_.add_inequality_constraint(quadratic_constraint<n>, &cons, 1e-8);
        }
    }
}

template<std::size_t n>
void base_optimizer<n>::results(std::vector<std::vector<double>>& dst)
{
    utils::copy_matvec2matvec(results_, dst);
}

template<std::size_t n>
std::vector<std::vector<double>> base_optimizer<n>::results()
{
    return results_;
}


template<std::size_t n>
void base_optimizer<n>::samples(std::vector<std::vector<double>>& dst)
{
    utils::copy_matvec2matvec(samples_, dst);
}

template<std::size_t n>
std::vector<std::vector<double>> base_optimizer<n>::samples()
{
    return samples_;
};  


template<std::size_t n>
void base_optimizer<n>::save_results(const std::string &name)
{
    utils::write_vec2file(results_, name);
};


template<std::size_t n>
void base_optimizer<n>::save_samples(const std::string &name)
{
    utils::write_vec2file(samples_, name);
};

// base_optimizer member functions finished
// biased_optimizer below

template<std::size_t n>
biased_optimizer<n>::biased_optimizer()
{
    std::cout <<"biased"  <<std::endl;
}

template<std::size_t n>
void biased_optimizer<n>::run(const int niter){
    std::vector<double> x(n);
    this->results_.resize(niter);
    this->samples_.resize(niter);
    for (int i=0; i<niter; ++i)
    {
        b_.x0 = this->uni_.sample();
        utils::copy_eig2vec(b_.x0, x);
        this->opt_.set_min_objective(biased_objective<n>, &b_);
        try
        {
            this->opt_.optimize(x, this->minf_);
            this->results_[i].resize(n);
            this->samples_[i].resize(n);
            utils::copy_vec2vec(x, this->results_[i]);
            utils::copy_eig2vec(b_.x0, this->samples_[i]);
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
slack_optimizer<n,m,l>::slack_optimizer(): base_optimizer<n+m+l+l>::base_optimizer(){};


template<std::size_t n, std::size_t m, std::size_t l>
slack_optimizer<n,m,l>::slack_optimizer(
            constraint_coeffs<n>& cons,
            const std::vector<double>& lb,
            const std::vector<double>& ub)
            : base_optimizer<n+m+l+l>::base_optimizer()
{
    this->set_bounds(lb, ub);
    this->add_constraints(cons);
};

template<std::size_t n, std::size_t m, std::size_t l>
slack_optimizer<n,m,l>::slack_optimizer(
            constraint_coeffs<n>& ineqcons,
            constraint_coeffs<n>& eqcons,
            const std::vector<double>& lb,
            const std::vector<double>& ub)
            : base_optimizer<n+m+l+l>::base_optimizer()
{
    this->set_bounds(lb, ub);
    this->add_constraints(ineqcons);
    this->add_constraints(eqcons);
};


template<std::size_t n, std::size_t m, std::size_t l>
void slack_optimizer<n,m,l>::set_bounds(const std::vector<double>& lb, const std::vector<double>& ub)
{
    std::vector<double> lb_new;
    std::vector<double> ub_new;
    utils::copy_vec2vec(lb, lb_new);
    utils::copy_vec2vec(ub, ub_new);
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
    base_optimizer<n+m+l+l>::set_bounds(lb_new, ub_new);
}

template<std::size_t n, std::size_t m, std::size_t l>
void slack_optimizer<n,m,l>::add_constraints(constraint_coeffs<n>& cons)
{
    assert (cons.type == "ineq" || cons.type == "eq");
    constraint_coeffs<n+m+l+l> cons_new;
    
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
        base_optimizer<n+m+l+l>::add_constraints(ineq_cons_ex_[ineq_count_-1]);
    }
    else 
    {
        cons_new.coeffs(n+m+eq_count_) = -1;
        constraint_coeffs<n+m+l+l> cons_new_eq_up;
        constraint_coeffs<n+m+l+l> cons_new_eq_low;
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
        base_optimizer<n+m+l+l>::add_constraints(eq_cons_ex_[eq_count_-1]);
        base_optimizer<n+m+l+l>::add_constraints(eq_cons_ex_up_[eq_count_-1]);
        base_optimizer<n+m+l+l>::add_constraints(eq_cons_ex_low_[eq_count_-1]);
    }
}



template<std::size_t n, std::size_t m, std::size_t l>
double slack_optimizer<n,m,l>::find_slack(const std::vector<double>& x, constraint_coeffs<n>& coeffs)
{
    double slack;
    std::vector<double> grad;
    if (coeffs.constype == "linear")
    {
        slack = linear_constraint<n>(x, grad, &coeffs);
    }
    else
    {
        slack = quadratic_constraint<n>(x, grad, &coeffs);
    }
    return slack;
}


template<std::size_t n, std::size_t m, std::size_t l>
double slack_optimizer<n,m,l>::find_slack(const std::vector<double>& x, constraint_coeffs<n>* coeffs)
{
    double slack;
    std::vector<double> grad;
    if (coeffs->constype == "linear")
    {
        slack = linear_constraint<n>(x, grad, coeffs);
    }
    else
    {
        slack = quadratic_constraint<n>(x, grad, coeffs);
    }
    return slack;
}


template<std::size_t n, std::size_t m, std::size_t l>
void slack_optimizer<n,m,l>::sample(std::vector<double> &x)
{
    // loop over ineq constraint coeffs find slack
    constraint_coeffs<n>*  coeffs_temp;
    double slack{0}, slack_temp{0};
    std::vector<double> x0;
    utils::copy_eig2vec(this->uni_.sample(),x);
    utils::copy_vec2vec(x, x0);
    x0.resize(n);
    //utils::copy_eig2vec(this->uni_.sample().head(n),x0);
    // initialize slack for inequalities
    for (int i=0; i<m; i++)
    {
        coeffs_temp = ineq_cons_[i];
        slack_temp = this->find_slack(x0 ,coeffs_temp);
        slack_temp = (slack_temp < 0) ? 0 : slack_temp;
        slack = (slack_temp > slack) ? slack_temp : slack;
        x[n+i] = slack;
            
    };
    slack = 0;
    for (int i=0; i<l; i++)
    {
        coeffs_temp = eq_cons_[i];
        slack = this->find_slack(x0 ,coeffs_temp);
        x[n+m+i] = slack;
        // choose t_i to be abs(slack_i)
        x[n+m+l+i] = (slack < 0) ? slack*-1 : slack;
    };
}


template<std::size_t n, std::size_t m, std::size_t l>
void slack_optimizer<n,m,l>::run(const int& niter)
{
    init_slack_data<n,m,l>(up_);
    std::vector<double> x(n+m+l+l);
    this->results_.resize(niter);
    this->samples_.resize(niter);
    for (int i=0; i<niter; ++i)
    {
        this->sample(x);
        this->samples_[i].resize(n+m+l+l);
        utils::copy_vec2vec(x, this->samples_[i]);
        this->opt_.set_min_objective(slack_objective<n,m,l>, &up_);
        try
        {
            this->opt_.optimize(x, this->minf_);
            this->results_[i].resize(n);
            utils::copy_vec2vec(x, this->results_[i]);
        }
        catch(std::exception &e)
        {
            std::cerr << e.what() << std::endl;
        }
    }
}

#endif



