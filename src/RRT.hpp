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
    bool valid{true};
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
        void addNode(node<n> q_new, int nearest_int);
        node<n> getNode(int id);
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
node<n> tree<n>::getNode(int id){
    return nodes_[id];
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
void tree<n>::addNode(node<n> q_new, int nearest_ind){
    node<n> q_near = nodes_[nearest_ind];
    addNode(q_new, q_near);
}

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


template< std::size_t n, std::size_t m>
class RRT : public BaseSampler<n,m> {
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
        void setRoot(node<n>);
        std::pair<node<n>,int> sampleNewNode();
        node<n> getNewNode(const node<n> &q_near, const node<n> &q_target, const double &alpha);
        virtual void run(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub);
        virtual void run(int n_iter);
        virtual void runOnTangent(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub);
        void run(int n_iter, Vector &lb, Vector &ub);
        void addNode(node<n>, int nearest_ind);
        bool tangentSampleInBounds(const Matrix<double,n-m,1> &u, TangentSpace<n,m> &tang);
    //    std::vector<std::vector<double>> results();
   //     void saveResults(std::string name);
     //   void saveSamples(std::string name);

    private:
        // step size
        bool use_tangent_{false};
        TargetProb<n,m> feasible_region_;
        double alpha_;
        UniformSampler<n,m> uni_;
        tree<n> tree_;
        double h_{1e-8}; // step width for numeric jacobian evaluation
};


template< std::size_t n, std::size_t m>
template<typename T>
void RRT<n,m>::setBounds(const T &lb, const T &ub)
{
    BaseSampler<n,m>::setBounds(lb, ub);
    uni_.setBounds(lb, ub);
}

template< std::size_t n, std::size_t m>
template<typename T>
RRT<n,m>::RRT(const T &lb, const T &ub, const double alpha, bool use_tangent){
    setBounds(lb,ub);
    setStepSize(alpha);
    setUseTangent(use_tangent);
}

template< std::size_t n, std::size_t m>
void RRT<n,m>::setUseTangent(bool use){
    use_tangent_ = use;
}

template< std::size_t n, std::size_t m>
void RRT<n,m>::setStepSize(const double &alpha){
    alpha_ = alpha;
}

template< std::size_t n, std::size_t m>
void RRT<n,m>::setH(double h){
    h_ = h;
}

template<std::size_t n, std::size_t m>
void RRT<n,m>::setRoot(node<n> root){
    tree_.setRoot(root);
}


/*template< std::size_t n, std::size_t m>
void RRT<n>::addConstraint(const ConstraintCoeffs<n> &con)
{
    feasible_region_.cons.push_back(con);
}*/

template<std::size_t n, std::size_t m>
std::pair<node<n>, int> RRT<n,m>::sampleNewNode(){
    node<n> target_node, new_node, nearest_node; 
    target_node.location = uni_.sample();
    nearest_node = tree_.findNearestNode(target_node);
    if (squaredDistNodes(nearest_node, target_node) > alpha_) {
        new_node = getNewNode(nearest_node, target_node, alpha_);
        if (!this->checkFeasible(new_node.location)) {
            new_node.valid = false;
        }
    } else {
        new_node.valid = false;
    }
    return std::make_pair(new_node,nearest_node.id);
}


template< std::size_t n, std::size_t m>
node<n> RRT<n,m>::getNewNode(const node<n> & q_near, const node<n> &q_target, const double &alpha)
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
template< std::size_t n, std::size_t m>
bool RRT<n>::checkFeasible(const Vector &q_new)
{
    return (feasible_region_(q_new) == 1) ? true: false;
};*/

template<std::size_t n, std::size_t m>
bool RRT<n,m>::tangentSampleInBounds(const Matrix<double,n-m,1> &u, TangentSpace<n,m> &tang){
    Matrix<double,n,1> lb(this->lb_), ub(this->ub_), u1;
    u1 = tang.toAmbient(u);
    return boundsCheck<n>(u1,lb,ub);
};

template<std::size_t n, std::size_t m>
void RRT<n,m>::addNode(node<n> new_node, int nearest_ind){
    tree_.addNode(new_node, nearest_ind);
    std::vector<double> new_node_vec(n);
    utils::copyEig2Vec(new_node.location, new_node_vec);
    this->results_.push_back(new_node_vec);
}


template< std::size_t n, std::size_t m>
void RRT<n,m>::run(int n_iter, Vector &lb, Vector &ub){
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
        nearest_node  = tree_.findNearestNode(target_node);
//        if (squaredDistNodes(nearest_node, target_node) > alpha_) {
            new_node = getNewNode(nearest_node, target_node, alpha_);
            if (this->checkFeasible(new_node.location) && boundsCheck<n>(new_node.location,lb,ub)){
                tree_.addNode(new_node, nearest_node);
                utils::copyEig2Vec(new_node.location, new_node_vec);
                this->results_.push_back(new_node_vec);
            }
  //      }
        this->samples_.push_back(utils::copyEig2Vec(target_node()));
    }
}

/*
template< std::size_t n, std::size_t m>
std::vector<std::vector<double>> RRT<n,m>::results(){
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
};*/

/*
template< std::size_t n, std::size_t m>
void RRT<n>::saveResults(std::string name){
    utils::writeVec2File(results(), name);
};


template< std::size_t n, std::size_t m>
void RRT<n>::saveSamples(std::string name){
    utils::writeVec2File(this->samples_, name);
};*/

template< std::size_t n, std::size_t m>
void RRT<n,m>::run(int n_iter){
    run(n_iter, this->lb_, this->ub_);
}


template< std::size_t n, std::size_t m>
void RRT<n,m>::runOnTangent(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub){
    //this->checkNumConstraints();
    // vector of functions that eval constraint
    // get jacobian in location x

    assert (this->cons_ptr_.size() == m);
    assert (m > 0);
    if (m==0) exit(1);
    std::vector<std::function<double(std::vector<double>&)>> funcs;
    ConstraintCoeffs<n> * c_ptr; 
    for (int i=0; i<this->cons_ptr_.size(); i++){
        c_ptr = this->cons_ptr_[i]; 
        if (c_ptr->type == "eq"){
            auto f = std::bind(evaluateConstraint<n>,_1, *c_ptr);
            funcs.push_back(f);
        }
    }
    Matrix<double,m,n> jac = numericJacobian<n,m>(seed,funcs,h_);
    TangentSpace<n,m> tang;

    Matrix<double,n,1> x0(seed.data());
    tang.findTangentSpace(x0,jac);
    //set lb and ub size to fit smaller rrt tang space
    std::vector<double> lb_ambient, ub_ambient;
    lb_ambient = utils::slice(lb,0,n-m-1);
    ub_ambient = utils::slice(ub,0,n-m-1);
    //RRT<n-m,0> rrt_tang(lb_ambient,ub_ambient,alpha_,false);
    //node<n-m> root_tangent, target_node, new_node;
    //root_tangent.location = Matrix<double,n-m,1>::Zero();
    //rrt_tang.setRoot(root_tangent);
    
    // set root in ambient space
    tree<n> tree_tangent_ambient; // tree capturing the tangent within the ambient space
    node<n> nearest_node, root_ambient, target_node, new_node_tangent, new_node_manifold;
    root_ambient.location = x0;
    this->setRoot(root_ambient);
    tree_tangent_ambient.setRoot(root_ambient);
    //node<n> sample_ambient_node, result_ambient_node;
    //node<n-m> sample_tangent_node, sample_tangent_node_nearest;
    int nearest_ind;
    std::vector<double> x_ambient(n);
    //std::pair<node<n-m>,int> sample_pair;

    UniformSampler<n-m,0> uni_local(lb_ambient, ub_ambient);

    for (int i=0; i<n_iter; i++){
        target_node.location = tang.toAmbient(uni_local.sample());
        nearest_node  = tree_tangent_ambient.findNearestNode(target_node);
        if (squaredDistNodes(nearest_node, target_node) > alpha_) {
            new_node_tangent = getNewNode(nearest_node, target_node, alpha_);
            if (boundsCheck<n>(new_node_tangent.location,this->lb_,this->ub_)){
                tree_tangent_ambient.addNode(new_node_tangent, nearest_node);
                utils::copyEig2Vec(new_node_tangent.location, x_ambient);
                this->samples_.push_back(x_ambient);
                x_ambient = this->optimize(x_ambient);
                new_node_manifold.location = Matrix<double,n,1>(x_ambient.data());
                nearest_node = tree_.getNode(nearest_node.id);
                tree_.addNode(new_node_manifold, nearest_node);
                this->results_.push_back(x_ambient);
            }
        this->samples_.push_back(utils::copyEig2Vec(target_node()));
        }
    }/*
    //    sample_pair = rrt_tang.sampleNewNode();
        sample_tangent_node = sample_pair.first;
        nearest_ind = sample_pair.second;
        sample_ambient_node.location = tang.toAmbient(sample_tangent_node.location);
        this->samples_.push_back(utils::copyEig2Vec(sample_ambient_node.location));
        if (sample_tangent_node.valid && tangentSampleInBounds(sample_tangent_node.location,tang)){
            rrt_tang.addNode(sample_tangent_node,nearest_ind);
            tree_tangent_ambient.addNode(sample_ambient_node, nearest_ind);
            utils::copyEig2Vec(sample_ambient_node.location, x_ambient); 
            this->samples_.push_back(x_ambient);
            x_ambient = this->optimize(x_ambient);
            result_ambient_node.location = Matrix<double,n,1>(x_ambient.data());    
            this->addNode(result_ambient_node, nearest_ind);
        }
   i*/ 
}

template< std::size_t n, std::size_t m>
void RRT<n,m>::run(int n_iter, std::vector<double> &seed, std::vector<double> &lb, std::vector<double> &ub){
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

  




#endif
