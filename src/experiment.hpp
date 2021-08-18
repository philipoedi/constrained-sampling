#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <memory>
#include <cmath>
#include "optimizer.hpp"
#include <cassert>
#include <iostream>
#include <vector>
#include "constraints.hpp"
#include "sampler.hpp"
#include "objectives.hpp"
#include <algorithm>
#include "RRT.hpp"
#include "kde.hpp"
#include "filter.hpp"
#include <filesystem>

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
        void setLocalWidths(const std::vector<double> &widths);
        void setGlobalWidths(const std::vector<double> &widths);
        void addConstraints(const ConstraintCoeffs<n> &con);
        void addConstraints(const std::vector<ConstraintCoeffs<n>> &cons);
        bool validOptimizer(std::string opt);
        bool validSampler(std::string samp);
        void setBandwidth(double bandwidth);
        void setBandwidthEstimator(std::string estimator);
        void setGridSpacing(double delta);
        void setSphere(double r);
        void setRBallWalk(double r);
        void setSave(bool save);
        void setFilter(double d_min);    
        void run();
        void setSuffix(std::string suf);
    

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
        double d_min_;
        bool global_use_tangent_{false};
        bool local_use_tangent_{false};
        // run parameters
        int global_n_iter_;
        int local_n_iter_;
        // space configurations
        std::vector<ConstraintCoeffs<n>> cons_;
        // output
        std::vector<std::vector<double>> local_results_;
        std::vector<std::vector<double>> global_results_;
        std::vector<std::vector<double>> global_samples_;
        std::vector<std::vector<double>> local_samples_;
        
        std::vector<double> widths_local_;
        std::vector<double> widths_global_;

        std::string suffix_{""};
        // kernel density estimator
        double delta_{1e-3};
        double sphere_radius_{0};
        double r_ballwalk_{0};
        double bandwidth_{0.1};
        std::string bandwidth_estimator_{"silverman"};
        bool use_bandwidth_estimator_{false};
        bool use_local_optimizer_{true};
        bool use_global_optimizer_{true};
        bool save_{false};
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
    use_local_optimizer_ = true;
    local_optimizer_ = local_optimizer;
}


template<std::size_t n, std::size_t m>
void Experiment<n,m>::setLocalSampler(std::string local_sampler){
    local_sampler_ = local_sampler;
}


template<std::size_t n, std::size_t m>
void Experiment<n,m>::setGlobalOptimizer(std::string global_optimizer){
    use_global_optimizer_ = true;
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
void Experiment<n,m>::setSave(bool save){
    save_ = save;
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
void Experiment<n,m>::setFilter(double d_min) {
    d_min_ = d_min;
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setSuffix(std::string suf){
    suffix_ = suf;
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
void Experiment<n,m>::setLocalWidths(const std::vector<double> &widths){
    widths_local_ = widths;
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setGlobalWidths(const std::vector<double> &widths){
    widths_global_ = widths;
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
void Experiment<n,m>::setGridSpacing(double delta){
    delta_ = delta;
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setBandwidth(double bandwidth){
    bandwidth_ = bandwidth;
    use_bandwidth_estimator_ = false;
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setBandwidthEstimator(std::string estimator){
    use_bandwidth_estimator_ = true;
    bandwidth_estimator_ = estimator;
}


template<std::size_t n, std::size_t m>
void Experiment<n,m>::setSphere(double r){
    sphere_radius_ = r;
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::setRBallWalk(double r){
    r_ballwalk_ = r;
}

template<std::size_t n, std::size_t m>
bool Experiment<n,m>::validOptimizer(std::string opt){
    if (opt == "rejection" || opt == "biased" || opt == "slack" || opt == "" || opt == "sphere" || opt == "circle" || opt == "line"){
        return true;
    }
    else {
        std::cout << "Choose any of the following optimizers: 'slack', 'biased' or if creating a reference dataset 'circle', 'sphere', 'line', 'rejection'" << std::endl;
        return false;
    }
}


template<std::size_t n, std::size_t m>
bool Experiment<n,m>::validSampler(std::string samp){
    if (samp == "RRT" || samp == "grid-walk" || samp == "uniform" || samp == "" || samp == "reference"){
        return true;
    }
    else {
        std::cout << "Choose any of the following samplers: 'uniform', 'grid-walk', 'RRT' or 'reference'" << std::endl;
        return false;
    }
}

template<std::size_t n, std::size_t m>
void Experiment<n,m>::run(){
    assert (validSampler(local_sampler_)); 
    assert (validSampler(global_sampler_)); 
    assert (validOptimizer(global_optimizer_)); 
    assert (validOptimizer(local_optimizer_)); 
    
    std::unique_ptr<BaseSampler<n,m>> global_sampler_ptr{nullptr};
    std::unique_ptr<BaseSampler<n,m>> local_sampler_ptr{nullptr};
    std::shared_ptr<BaseOptimizer<n>> global_optimizer_ptr{nullptr};
    std::shared_ptr<BaseOptimizer<n>> local_optimizer_ptr{nullptr};
 
    std::string name, folder;
    name = utils::getDateTimeString();
    name = name + "_" + global_sampler_ + "_" + global_optimizer_+ "_" + local_sampler_ +"_" + local_optimizer_;

    // create folder for results
    folder = name + suffix_;
    std::filesystem::create_directory(folder);
    name = folder + "/" + name;
 
    int last_n_iter{local_n_iter_};
    
    int total_n_iter = (local_n_iter_ == 0) ? global_n_iter_ : global_n_iter_ * local_n_iter_;

    if (global_sampler_  == "reference"){ 
        if (global_optimizer_  == "sphere"){
            std::vector<double> lb_sphere{0,-sphere_radius_};
            std::vector<double> ub_sphere{2*M_PI,sphere_radius_};
            SphereSampler ssamp(lb_sphere, ub_sphere);
            ssamp.run(total_n_iter);
            ssamp.saveSamples(name+"_global");
            ssamp.saveResults(name+"_global");
            ssamp.saveResults(name+"_local_0");
            ssamp.saveResults(name+"_local_0");
            local_samples_ = ssamp.results();
        } else if ( global_optimizer_ == "circle") {
            std::vector<double> lb_circle{0};
            std::vector<double> ub_circle{2*M_PI};
            CircleSampler csamp(lb_circle, ub_circle);
            csamp.setRadius(sphere_radius_);
            csamp.run(total_n_iter);
            csamp.saveSamples(name+"_global");
            csamp.saveResults(name+"_global");
            csamp.saveResults(name+"_local_0");
            csamp.saveResults(name+"_local_0");
            local_samples_ = csamp.results();
    
        } else if (global_optimizer_  == "line") {
            LineSampler lsamp(lb_global_, ub_global_);
            ConstraintCoeffs<2> line_con;
            line_con.coeffs << cons_[0].coeffs(0), cons_[0].coeffs(1);
            line_con.cons = cons_[0].cons;
            std::vector<ConstraintCoeffs<2>> line_cons_vec;
            line_cons_vec.push_back(line_con);
            lsamp.addConstraints(line_cons_vec);
            lsamp.run(total_n_iter);        
            lsamp.saveSamples(name+"_global");
            lsamp.saveResults(name+"_global");
            lsamp.saveResults(name+"_local_0");
            lsamp.saveResults(name+"_local_0");
            local_samples_ = lsamp.results();
        } else {
            UniformSampler<n,m> usamp(lb_global_,ub_global_);
            usamp.addConstraints(cons_);
            usamp.run(total_n_iter);
            usamp.saveSamples(name+"_global");
            usamp.saveResults(name+"_global");
            usamp.saveResults(name+"_local_0");
            usamp.saveResults(name+"_local_0");
            local_samples_ = usamp.results();
        }
    } else {
        if (global_sampler_ == "uniform"){
            global_sampler_ptr = std::make_unique<UniformSampler<n,m>>(lb_global_, ub_global_);
        } else if (global_sampler_ == "RRT") {
            global_sampler_ptr = std::make_unique<RRT<n,m>>(lb_global_, ub_global_, global_alpha_, global_use_tangent_);
        } else if (global_sampler_ == "grid-walk") {
            global_sampler_ptr = std::make_unique<GridWalk<n,m>>(lb_global_, ub_global_, widths_global_, global_use_tangent_);
        } else {
            global_sampler_ptr = std::make_unique<GridWalk<n,m>>(lb_global_, ub_global_, widths_global_, global_use_tangent_, true);
        }
       /*
        if (global_optimizer_ == "biased") {
            //global_optimizer_ptr = std::make_shared<BiasedOptimizer<n>>(lb_global_, ub_global_);
            global_sampler_ptr->setOptimizer("biased",);
        }*/

       // global sampler
        global_sampler_ptr->addConstraints(cons_);
        if (use_global_optimizer_) global_sampler_ptr->setOptimizer(global_optimizer_, lb_global_, ub_global_);
        //if (global_optimizer_) global_sampler_ptr->setOptimizer(global_optimizer_ptr);

        // global results
        global_sampler_ptr->run(global_n_iter_);
        global_samples_ = global_sampler_ptr->samples();
        global_results_ = global_sampler_ptr->results();
        if (save_){
            global_sampler_ptr->saveResults(name+"_global");
            global_sampler_ptr->saveSamples(name+"_global");
            global_sampler_ptr->saveNumIterations(name+"_global");
        }

        std::vector<bool> accepted_results;
        if (d_min_ > 0) {
            std::cout << "global_results_.size()" << global_results_.size() << std::endl;
            int num_accepted{0};
            greedyNodeRemoval(global_results_, d_min_, accepted_results, num_accepted);
            local_n_iter_ = std::floor(total_n_iter / num_accepted);
            std::cout << "total_n_iter: " << total_n_iter << std::endl;
            std::cout << "local_n_iter: " << local_n_iter_ << std::endl;
            last_n_iter = local_n_iter_ + total_n_iter - (local_n_iter_ * global_results_.size());
            std::cout << "global_results_.size()" << global_results_.size() << std::endl;
        }

        bool local{true};

        if (local_sampler_ == "uniform") {
            local_sampler_ptr = std::make_unique<UniformSampler<n,m>>(lb_global_, ub_global_);
        } else if (local_sampler_ == "RRT"){
            local_sampler_ptr = std::make_unique<RRT<n,m>>(lb_global_, ub_global_, local_alpha_, local_use_tangent_);
        } else if (local_sampler_ == "grid-walk") {
            local_sampler_ptr = std::make_unique<GridWalk<n,m>>(lb_global_, ub_global_, widths_local_, local_use_tangent_);
        } else if (local_sampler_ == "ball-walk"){
            // case ballwalk
            local_sampler_ptr = std::make_unique<GridWalk<n,m>>(lb_global_, ub_global_, widths_local_, local_use_tangent_, r_ballwalk_);
        } else {
            local_samples_ = global_results_;
            local = false;
        }
      /* 
        if (local_optimizer_ == "biased"){
            local_optimizer_ptr = std::make_shared<BiasedOptimizer<n>>(lb_global_, ub_global_);
        }*/

        
        if (local) {
            // local sampler
            local_sampler_ptr->addConstraints(cons_);
            if (use_local_optimizer_) local_sampler_ptr->setOptimizer(local_optimizer_, lb_global_, ub_global_);

            // local results 
            for (int i=0; i<global_results_.size() ;i++){
                std::cout << i << std::endl;
                if (d_min_ > 0 && !accepted_results[i]) continue; 
                if (i-1 == global_results_.size()){ 
                    local_n_iter_ = last_n_iter;
                }
                if (local_use_tangent_) {
                   local_sampler_ptr->runOnTangent(local_n_iter_, global_results_[i], lb_local_, ub_local_);
                } else {
                   local_sampler_ptr->run(local_n_iter_, global_results_[i], lb_local_, ub_local_);
                }
                if (save_) {
                    local_sampler_ptr->saveSamples(name+"_local_"+std::to_string(i));
                    local_sampler_ptr->saveResults(name+"_local_"+std::to_string(i));
                    local_sampler_ptr->saveNumIterations(name+"_local_"+std::to_string(i));
                }
                utils::appendVec2Vec(local_sampler_ptr->results(),local_samples_);
                local_sampler_ptr->reset();
            }
        }
        // probability density estimation of local samples

    }
    KernelEstimator<n,n> kdest;
    kdest.fit(local_samples_);
    if (use_bandwidth_estimator_) {
        kdest.find_optimal_bandwidth(bandwidth_estimator_);
    } else {
        kdest.setBandwidth(bandwidth_);
    }
    if (sphere_radius_ > 0 && global_optimizer_ != "circle") {
        kdest.setSphere(sphere_radius_);
    }
    kdest.predict(lb_global_ ,ub_global_ ,delta_);
    if (save_) {
        kdest.savePdes(name+"_pdes");
    }
    std::vector<double> probs;
    probs = kdest.leaveOneOutEstimation();
    if (save_) {
        utils::writeVec2File(probs,name+"_probs");
    }
    //utils::shuffleVector(local_samples_);
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
