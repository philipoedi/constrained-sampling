#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "optimizer.hpp"
#include <cassert>
#include <iostream>
#include <vector>
#include "constraints.hpp"
#include "sampler.hpp"
#include "objectives.hpp"
#include <algorithm>
#include "RRT.hpp"

template<std::size_t n, std::size_t m>
class Experiment {

    public:
        Experiment();
        Experiment(
            std::string global_optimizer,
            std::string global_sampler,
            std::string local_optimizer,
            std::string local_sampler);
        void setGlobalSampler(std::string global_sampler);
        void setLocalSampler(std::string local_sampler);
        void setGlobalOptimizer(std::string global_optimizer);
        void setLocalOptimizer(std::string local_optimizer);
        void setGlobalAlpha(double global_alpha);
        void setLocalAlpha(double local_alpha);
        void setGlobalBounds(const std::vector<double> &lb, const std::vector<double> &ub);
        void setLocalBounds(const std::vector<double> &lb, const std::vector<double> &ub);
        void setLocalLowerBounds(const std::vector<double> &lb);
        void setGlobalLowerBounds(const std::vector<double> &lb);
        void setLocalUpperBounds(const std::vector<double> &ub);
        void setGlobalUpperBounds(const std::vector<double> &ub);
        void setLocalNumIter(int local_num_iter);
        void setGlobalNumIter(int global_num_iter);
        void setLocalUseTangent(bool use);
        void setGlobalUseTangent(bool use);
        void addConstraints(const ConstraintCoeffs<n> &con);
        void addConstraints(const std::vector<ConstraintCoeffs<n>> &cons);
        bool validOptimizer(std::string opt);
        bool validSampler(std::string samp);
        void run();

    private:
        std::string global_sampler_; // "uniform","rrt"
        std::string local_sampler_; // "uniform","normal","RRT","metropolis_hastings"
        std::string global_optimizer_; // "slack","biased", "none"
        std::string local_optimizer_; // "slack", "biased", "tangent_normal"
        std::vector<double> lb_global_;
        std::vector<double> ub_global_;
        std::vector<double> lb_local_;
        std::vector<double> ub_local_;
        // RRT related parameters
        double global_alpha_;
        double local_alpha_;
        bool global_use_tangent_;
        bool local_use_tangent_;
        // run parameters
        int global_n_iter_;
        int local_n_iter_;
        // space configurations
        std::vector<ConstraintCoeffs<n>> cons_;
        // output
        std::vector<std::vector<double>> local_results_;
        std::vector<std::vector<double>> global_results_;
        std::vector<std::vector<double>> global_samples_;
        std::vector<std::vector<double>> local_samples;

};

template<std::size_t n, std::size_t m>
Experiment<n,m>::Experiment(){

};


template<std::size_t n, std::size_t m>
Experiment<n,m>::Experiment(
    std::string global_optimizer,
    std::string global_sampler,
    std::string local_optimizer,
    std::string local_sampler):
    global_sampler_(global_sampler),
    global_optimizer_(global_optimizer),
    local_sampler_(local_sampler),
    local_optimizer_(local_optimizer) {

}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setLocalOptimizer(std::string local_optimizer){
    local_optimizer_ = local_optimizer;
}


template<std::size_t n, std::size_t m>
void Experiment<n,m>::setLocalSampler(std::string local_sampler){
    local_sampler_ = local_sampler;
}


template<std::size_t n, std::size_t m>
void Experiment<n,m>::setGlobalOptimizer(std::string global_optimizer){
    global_optimizer_ = global_optimizer;
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setGlobalSampler(std::string global_sampler){
    global_sampler_ = global_sampler;
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setLocalLowerBounds(const std::vector<double> &lb){
    lb_local_ = lb;
}


template<std::size_t n, std::size_t m>
void Experiment<n,m>::setGlobalLowerBounds(const std::vector<double> &lb){
    lb_global_= lb;
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setLocalUpperBounds(const std::vector<double> &ub){
    ub_local_ = ub;
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setGlobalUpperBounds(const std::vector<double> &ub){
    ub_global_ = ub;
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setGlobalBounds(const std::vector<double> &lb, const std::vector<double> &ub){
    setGlobalUpperBounds(ub);
    setGlobalLowerBounds(lb);
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setLocalBounds(const std::vector<double> &lb, const std::vector<double> &ub){
    setLocalUpperBounds(ub);
    setLocalLowerBounds(lb);
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setGlobalAlpha(double alpha){
    global_alpha_ = alpha;
}

template<std::size_t n,std::size_t m>
void Experiment<n,m>::setLocalAlpha(double alpha){
    local_alpha_ = alpha;
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setLocalUseTangent(bool use){
    local_use_tangent_ = use;
}


template<std::size_t n, std::size_t m>
void Experiment<n,m>::setGlobalUseTangent(bool use){
    global_use_tangent_ = use;
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setGlobalNumIter(int num_iter){
    global_n_iter_ = num_iter;
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setLocalNumIter(int num_iter){
   local_n_iter_ = num_iter;
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::addConstraints(const ConstraintCoeffs<n> &con){
    cons_.push_back(con);
}


template<std::size_t n, std::size_t m>
void Experiment<n,m>::addConstraints(const std::vector<ConstraintCoeffs<n>> &cons){
    std::copy(cons.begin(),cons.end(),std::back_inserter(cons_));
}

template<std::size_t n, std::size_t m>
bool Experiment<n,m>::validOptimizer(std::string opt){
    if (opt == "biased" || opt == "slack" || "none"){
        return true;
    }
    else {
        std::cout << "Choose any of the following optimizers: 'slack', 'biased'" << std::endl;
        return false;
    }
}


template<std::size_t n, std::size_t m>
bool Experiment<n,m>::validSampler(std::string samp){
    if (samp == "RRT" || samp == "metropolis_hastings" || samp == "uniform"){
        return true;
    }
    else {
        std::cout << "Choose any of the following samplers: 'uniform', 'metropolis_hastings', 'RRT'" << std::endl;
        return false;
    }
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::run(){
    assert (validSampler(local_sampler_)); 
    assert (validSampler(global_sampler_)); 
    assert (validOptimizer(global_optimizer_)); 
    assert (validOptimizer(local_optimizer_)); 
    
    BaseSampler<n> * global_sampler_ptr{nullptr};
    BaseSampler<n> * local_sampler_ptr{nullptr};
    BaseOptimizer<n> * global_optimizer_ptr{nullptr};
    BaseOptimizer<n> * local_optimizer_ptr{nullptr};
 

    if (global_sampler_ == "uniform"){
        global_sampler_ptr = new UniformSampler<n>(lb_global_, ub_global_);
    } else if (global_sampler_ == "RRT") {
        global_sampler_ptr = new RRT<n>(lb_global_, ub_global_, global_alpha_, global_use_tangent_);
    }
   
    if (global_optimizer_ == "biased") {
        global_optimizer_ptr = new BiasedOptimizer<n>(lb_global_, ub_global_);
    }

    std::string name;
    name = utils::getDateTimeString();
    name = name +"_" + global_sampler_ + "_" + global_optimizer_+ "_" + local_sampler_ +"_" + local_optimizer_;

    // global sampler
    global_sampler_ptr->addConstraints(cons_);
    global_sampler_ptr->setOptimizer(global_optimizer_ptr);

    // global results
    global_sampler_ptr->run(global_n_iter_);
    global_samples_ = global_sampler_ptr->samples();
    global_results_ = global_sampler_ptr->results();
    global_sampler_ptr->saveResults(name+"_global");
    global_sampler_ptr->saveSamples(name+"_global");


    if (local_sampler_ == "uniform") {
        local_sampler_ptr = new UniformSampler<n>(lb_global_, ub_global_);
    } else if (local_sampler_ == "RRT"){
        local_sampler_ptr = new RRT<n>(lb_global_, ub_global_, local_alpha_, local_use_tangent_);
    }
   
    if (local_optimizer_ == "biased"){
        local_optimizer_ptr = new BiasedOptimizer<n>(lb_global_, ub_global_);
    }

    
    // local sampler
    local_sampler_ptr->addConstraints(cons_);
    local_sampler_ptr->setOptimizer(local_optimizer_ptr);

    // local results 
    for (int i=0; i<global_results_.size() ;i++){
        local_sampler_ptr->run(local_n_iter_, global_results_[i], lb_local_, ub_local_);
        local_sampler_ptr->saveSamples(name+"_local_"+std::to_string(i));
        local_sampler_ptr->saveResults(name+"_local_"+std::to_string(i));
        local_sampler_ptr->reset();
    }
    
}




/*
Experiment::run(){
    // global sampler

    

    switch(global_sampler) {
        case "uniform":
           :wq

        case "RRT":
        
        default:
            std::cout << "Wrong global sampling method chosen. Use 'uniform' or 'RRT' 
    }

}

std::fnuction< >local_sample;
std::function<> global_optmier
std::function<> local_optimzer;
std::function<> global_sample;

if ()
    local_sample = 
elif
    local_sample = 
elif

global_samples
global_results
local_resulst
local_samples

for (niter_global)
    init localsampler // set seed in mh, set root in rrt;
    localsampler(n_iter_local);
    local_results.append()
    local_samples.append()
*/

#endif
