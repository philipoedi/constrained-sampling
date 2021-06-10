#ifndef RRT_H
#define RRT_H

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
    Matrix<double,n,1> operator()(){ return location;};
};

template<std::size_t n>
Matrix<double,n,1> squaredDistNodes(const node<n> &n1, const node<n> &n2)
{
    return utils::squaredDist(n1(), n2());
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
void tree<n>::setRoot(node<n> root){
    root.id = 0;
    num_nodes_ += 1;
    root_ = root;
};

template<std::size_t n>
void tree<n>::addNode(node<n> q_new, const node<n> &q_near)
{
    q_new.id = num_nodes_;
    tree_.push_back(std::make_pair<>);
    nodes_.push_back(q_new); 
    num_nodes_+= 1; 
};

template<std::size_t n>
node<n> tree<n>::findNearestNode(const node<n> &q_new)
{
    std::vector<double> distances;
    int min_index;
    std::transform(nodes_.begin(),nodes_.end(), std::back_inserter(distances), std::bind(squaredDistNodes, q_new, _1));
    min_index = std::min_element(distances.begin(),distances.end()) - distances.begin();
    return nodes_[min_index]; 
};

template<std::size_t n>
void tree<n>::save(std::string name)
{
    std::vector<std::vector<double>> nodes_vec;
    std::transform(nodes_.begin(),nodes_.end(), std::back_inserter(nodes_vec.begin(), utils::copyEig2Vec);
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
        template<typename T> void setBounds(const T &lb, const T &ub);
        void setStepSize(const double &alpha);  
        node<n> getNewNode(const node<n> &q_near, const node<n> &q_target);
        bool checkQFree(const node<n> &q_new);
        void explore(int n_iter);
        void saveResults(std::string name);
        void saveSamples(std::string name);
        std::vector<std::vector<double>> results();

    private:
        // step size
        std::vector<std::vector<double>> samples_;
        double alpha_;
        UniformSampler<n> uni_;
        tree<n> tree_;
};


template<std::size_t n>
node<n> RRT<n>::getNewNode(const node<n> & q_near, const node<n> &q_target)
{
    node<n> q_new;
    Vector p_new, p_target, p_near;
    p_near = q_near();
    p_target = q_target();
    q_new.location = p_near + (alpha_ /  (p_target - p_near).squaredNorm().sqrt()) * (p_target - p_near);
    return q_new;
};

template<std::size_t n>
void RRT<n>::explore(int n_iter){
    node<n> target_node, nearest_node, new_node;
    bool found_root{false};
    while (!found_root && n_iter > 0){
        target_node.location = uni_.sample();
        if (checkFeasible(target_node())){
           tree_.setRoot(target_node);
           found_root = true; 
        }
        n_iter -= 1;
        samples_.push_back(utils::copyEig2Vec(target_node()));
    }
    for (int i=0; i< n_iter; i++)
    {
        target_node.location = uni_.sample();
        nearest_node  = tree_.findNearestNode(target_node);
        new_node = getNewNode(nearest_node);
        if (checkFeasible(new_node)){
            tree_.addNode(new_node);
        }
        samples_.push_back(utils::copyEig2Vec(target_node()));
    }
}

template<std::size_t n>
std::vector<std::vector<double>> RRT<n>::results(){
    std::vector<std::vector<double>> results(tree_.size(),n);
    std::vector<node<n>> nodes = tree_.getNodes();
    std::transform(nodes.begin(), nodes.end(), results.begin(), utils::copyEig2Vec);
    return results;
};

template<std::size_t n>
void RRT<n>::saveResults(std::string name){
    utils::writeVec2File(results(), name);
};


template<std::size_t n>
void RRT<n>::saveSamples(std::string name){
    utils::writeVec2File(samples_, name);
};




#endif
