#ifndef RRT_H
#define RRT_H

#include <functional>
#include "constraints.hpp"
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "sampler.hpp"
#include <utility>
#include "utils.hpp"


using namespace Eigen;
using namespace std::placeholders;


template<std::size_t n>
struct node {
    int id;
    Matrix<double,n,1> location;
    Matrix<double,n,1> operator()() const { return location;};
    Matrix<double,n,1> operator()() { return location;};
};

template<std::size_t n>
double squaredDistNodes(const node<n> &n1, const node<n> &n2)
{
    return utils::squaredDist<n>(n1(), n2());
};

template<std::size_t n>
std::vector<int> pair2Vec(const std::pair<node<n>,node<n>> p){
    std::vector<int> data;
    data.push_back(p.first());
    data.push_back(p.second());
    return data;  
};


template<std::size_t n>
class tree {
    
    typedef Matrix<double,n,1> Vector;
    enum {NeedsToAlign = (sizeof(Vector)%16) == 0};
    
    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign);
        tree(){};
        tree(const node<n> root);
        void setRoot(node<n> root);
        void addNode(node<n> q_new, const node<n> &q_near);
        node<n> findNearestNode(const node<n> &q_new);
        void save(std::string name);
        int size();
        std::vector<node<n>> getNodes();
        bool hasRoot();

    private:
        
        int num_nodes_{0}; 
        node<n> root_;
        std::vector<node<n>> nodes_;
        std::vector<std::pair<node<n>,node<n>>> tree_;

};

template<std::size_t n>
tree<n>::tree(const node<n> root){
    setRoot(root);
};

template<std::size_t n>
bool tree<n>::hasRoot(){
    return (num_nodes_ > 0) ? true : false;
}

template<std::size_t n>
int tree<n>::size() {return num_nodes_;};

template<std::size_t n>
std::vector<node<n>> tree<n>::getNodes() {return nodes_;};

template<std::size_t n>
void tree<n>::setRoot(node<n> root){
    root.id = 0;
    num_nodes_ += 1;
    root_ = root;
    nodes_.push_back(root);
};

template<std::size_t n>
void tree<n>::addNode(node<n> q_new, const node<n> &q_near)
{
    q_new.id = num_nodes_;
    std::pair<node<n>,node<n>> vertice;
    vertice = std::make_pair(q_near, q_new);
    tree_.push_back(vertice);
    nodes_.push_back(q_new); 
    num_nodes_+= 1; 
};

template<std::size_t n>
node<n> tree<n>::findNearestNode(const node<n> &q_new)
{
    std::vector<double> distances;
    distances.clear();
    int min_index;
    std::transform(nodes_.begin(),nodes_.end(), std::back_inserter(distances), std::bind(squaredDistNodes<n>, _1, q_new));;
    min_index = std::min_element(distances.begin(),distances.end()) - distances.begin();
    return nodes_[min_index]; 
};

template<std::size_t n>
void tree<n>::save(std::string name)
{
    std::vector<std::vector<double>> nodes_vec;
    typename std::vector<node<n>>::iterrator it;
    for (it = nodes_.begin(); it != nodes_.end(); ++it){
        nodes_vec.push_back(*it());
    };
    //std::transform(nodes_.begin(),nodes_.end(), std::back_inserter(nodes_vec.begin()), utils::copyEig2Vec);
    utils::writeVec2File(nodes_vec, name+"_nodes");
    std::vector<std::vector<int>> tree_vec;
    std::transform(tree_.begin(), tree_.end(), std::back_inserter(tree_vec), pair2Vec);
    utils::writeVec2File(tree_vec, name+"_tree");
};


template<std::size_t n>
class RRT : public BaseSampler<n> {
    typedef Matrix<double,n,1> Vector;
    enum {NeedsToAlign = (sizeof(Vector)%16) == 0};

    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign);
        RRT(){};
        template<typename T> RRT(const T &lb, const T &ub, const double alpha, bool use_tangent);
        template<typename T> void setBounds(const T &lb, const T &ub);
        void setUseTangent(bool use_tangent);
        void setStepSize(const double &alpha);  
        void setH(double h);
        node<n> getNewNode(const node<n> &q_near, const node<n> &q_target, const double &alpha);
        virtual void run(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub);
        virtual void run(int n_iter);
        void run(int n_iter, Vector &lb, Vector &ub);
        void runOnTangent(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub);
        std::vector<std::vector<double>> results();
   //     void saveResults(std::string name);
     //   void saveSamples(std::string name);

    private:
        // step size
        bool use_tangent_{false};
        TargetProb<n> feasible_region_;
        double alpha_;
        UniformSampler<n> uni_;
        tree<n> tree_;
        double h_{1e-8}; // step width for numeric jacobian evaluation
};


template<std::size_t n>
template<typename T>
void RRT<n>::setBounds(const T &lb, const T &ub)
{
    BaseSampler<n>::setBounds(lb, ub);
    uni_.setBounds(lb, ub);
}

template<std::size_t n>
template<typename T>
RRT<n>::RRT(const T &lb, const T &ub, const double alpha, bool use_tangent){
    setBounds(lb,ub);
    setStepSize(alpha);
    setUseTangent(use_tangent);
}

template<std::size_t n>
void RRT<n>::setUseTangent(bool use){
    use_tangent_ = use;
}

template<std::size_t n>
void RRT<n>::setStepSize(const double &alpha){
    alpha_ = alpha;
}

template<std::size_t n>
void RRT<n>::setH(double h){
    h_ = h;
}

/*template<std::size_t n>
void RRT<n>::addConstraint(const ConstraintCoeffs<n> &con)
{
    feasible_region_.cons.push_back(con);
}*/


template<std::size_t n>
node<n> RRT<n>::getNewNode(const node<n> & q_near, const node<n> &q_target, const double &alpha)
{
    node<n> q_new;
    Vector p_new, p_target, p_near;
    p_near = q_near();
    p_target = q_target();
    //std::cout << "target location" << q_target.location << std::endl;
    //std::cout << "nearest location" << q_near.location << std::endl;
    q_new.location = p_near + std::sqrt(alpha /  (p_target - p_near).squaredNorm()) * (p_target - p_near);
    //std::cout << "new location" << q_new.location << std::endl;
    return q_new;
};
/*
template<std::size_t n>
bool RRT<n>::checkFeasible(const Vector &q_new)
{
    return (feasible_region_(q_new) == 1) ? true: false;
};*/


template<std::size_t n>
void RRT<n>::run(int n_iter, Vector &lb, Vector &ub){
    node<n> target_node, nearest_node, new_node;
    std::vector<double> new_node_vec(n);
    while (!tree_.hasRoot()  && n_iter > 0){
        target_node.location = uni_.sample();
        if (this->checkFeasible(target_node.location) && boundsCheck<n>(target_node.location, lb, ub)){
           tree_.setRoot(target_node);
           std::cout << "Feasible starting node found" <<std::endl;
        }
        n_iter -= 1;
        this->samples_.push_back(utils::copyEig2Vec(target_node()));
    }
    for (int i=0; i< n_iter; i++)
    {
        target_node.location = uni_.sample();
        std::cout << "target node:" << target_node.location << std::endl;
        nearest_node  = tree_.findNearestNode(target_node);
        if (squaredDistNodes(nearest_node, target_node) > alpha_) {
            new_node = getNewNode(nearest_node, target_node, alpha_);
            if (this->checkFeasible(new_node.location) && boundsCheck<n>(new_node.location,lb,ub)){
                tree_.addNode(new_node, nearest_node);
                utils::copyEig2Vec(new_node.location, new_node_vec);
                this->results_.push_back(new_node_vec);
            }
        }
        this->samples_.push_back(utils::copyEig2Vec(target_node()));
    }
}

template<std::size_t n>
std::vector<std::vector<double>> RRT<n>::results(){
    std::vector<std::vector<double>> results;
    std::vector<node<n>> nodes = tree_.getNodes();
    typename std::vector<node<n>>::iterator it;
    Vector data;
    node<n> temp;
    for (it=nodes.begin(); it!=nodes.end(); ++it)
    {
        temp = *it;
        results.push_back(utils::copyEig2Vec(temp()));
    };
    return results;
};

/*
template<std::size_t n>
void RRT<n>::saveResults(std::string name){
    utils::writeVec2File(results(), name);
};


template<std::size_t n>
void RRT<n>::saveSamples(std::string name){
    utils::writeVec2File(this->samples_, name);
};*/

template<std::size_t n>
void RRT<n>::run(int n_iter){
    run(n_iter, this->lb_, this->ub_);
}

template<std::size_t n>
void RRT<n>::run(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub){
    Vector lb_old, ub_old, seed_eig(seed.data()), lb_eig(lb.data()), ub_eig(ub.data());
    lb_old = this->lb_;
    ub_old = this->ub_;
    setBounds(seed_eig + lb_eig, seed_eig + ub_eig);
    node<n> root;
    root.location = seed_eig;
    tree_.setRoot(root);
    run(n_iter, lb_old, ub_old);
    setBounds(lb_old, ub_old);
}


template<std::size_t n>
void RRT<n>::runOnTangent(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub){
    // vector of functions that eval constraint
    // get jacobian in location x
    std::vector<std::function<double(std::vector<double>&)>> funcs;
    for (int i=0; i<this->cons_ptr_.size(); i++){
        ConstraintCoeffs<n> * c_ptr = this->cons_ptr_[i]; 
        if (c_ptr->type == "eq"){
            auto f = std::bind(evaluateConstraint<n>,_1, *c_ptr);
            funcs.push_back(f);
        }
    }
    Matrix<double,funcs.size(),n> jac = numericJacobian<n,m>(seed,funcs,h_);

}
    //numericJacobian();
    // declare tang space
    // find tangspace
    
    // new rrt obj for tang space
    // create tree on tang space
    // check for each node if in ambient space within bounds
    // to ambient -> samples
    // all samples -> project using optimizer all samples -> project using optimizer::






#endif
