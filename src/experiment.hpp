#include <iostream>


class Experiment {

    public:
        setBounds;
        addConstraints;
        setGlobalSampler();
        setLocalSampler();
        setGlobalOptimizer();
        setLocalOptimizer();
        setGlobalAlpha();
        setLocalAlpha();
        setLocalBounds();
        setLocalNumIter();
        setGlobalNumIter();
        init();
        run();

    private:
        std::string global_sampler_; // "uniform","rrt"
        std::string local_sampler_; // "uniform","normal","RRT","metropolis_hastings"
        std::string global_optimizer_; // "slack","biased", "none"
        std::strnig local_optimizer_; // "slack", "biased", "tangent_normal"
        // global sampler
        // global optimizer
        // local sampler
        // local optimizer
        std::vector<ConstraintCoeffs> cons_;
        std::vector<double> lb_;
        std::vector<double> ub_;
        std::vector<double> lb_local_;
        std::vector<double> ub_local_;
        double global_alpha_;
        double local_alpha_;
        // other parameters
        std::vector<double> results_;
        int n_iter_local_;
        int n_iter_global_;
};



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


