#include "filter.hpp"



UpperTriMat::UpperTriMat(std::size_t n){
        resize(n);
}

std::size_t UpperTriMat::size(){
    return data_.size() + 1;
}

void UpperTriMat::resize(std::size_t n){
    data_.resize(n-1);
    for (int i=0; i<n-1; i++){
        data_[i].resize(i+1);
    }
}

double UpperTriMat::operator()(std::size_t i, std::size_t j){
    if (i > j) {
        return data_[i-1][j];
    } else if (i < j) {
        return data_[j-1][i];
    } else {
        return 0;
    }
}

std::vector<double> UpperTriMat::operator()(std::size_t i){
    return data_[i];
}


void UpperTriMat::addData(std::size_t i, std::size_t j, double x){
    if (i > j){
        data_[i-1][j] = x;
    } else {
        data_[j-1][i] = x;
    }
}


double distance(const std::vector<double> &a, const std::vector<double> &b){
    assert (a.size() == b.size());
    Map<const VectorXd> a_eig(a.data(), a.size());
    Map<const VectorXd> b_eig(b.data(), b.size());
    return (a_eig - b_eig).norm();    
}

UpperTriMat getDistanceMatrix(const std::vector<std::vector<double>> &samples){
    UpperTriMat d;
    double d_ij;
    std::size_t n = samples.size();
    d.resize(n);
    for (std::size_t i=0; i<n; i++){
        for (std::size_t j=0; j<n; j++){
            if (i>j){
                d_ij = distance(samples[i],samples[j]);
                d.addData(i,j,d_ij);
            }
        }
    }

    return d;
}

std::vector<std::pair<std::size_t,std::size_t>> filterVertices(UpperTriMat distances, double d_min){
    std::vector<double> row;
    std::vector<std::pair<std::size_t, std::size_t>> vertices_to_remove;
    for (std::size_t i=0; i<distances.size()-1; i++){
        row.clear();
        row = distances(i);
        for (std::size_t j=0; j<row.size(); j++){
            if (row[j] < d_min) {
                vertices_to_remove.push_back(std::make_pair(i+1,j));
            }
        }
    }
    return vertices_to_remove;
}


std::set<std::size_t> getUniqueNodes(std::vector<std::pair<std::size_t,std::size_t>> &vertices){
    std::set<std::size_t> unique_nodes;
    std::vector<std::pair<std::size_t,std::size_t>>::iterator it;
    std::pair<std::size_t,std::size_t> v;
    for (it=vertices.begin(); it!=vertices.end(); it++){
        v = *it;
        unique_nodes.insert(v.first);
        unique_nodes.insert(v.second);
    }
    return unique_nodes;
}

std::map<std::size_t,double> getNodeTotalDistanceMap(std::set<std::size_t> &nodes, UpperTriMat &distances){
    std::set<std::size_t>::iterator it;
    std::size_t n_i;
    double total_distance{0};
    std::map<std::size_t,double> node_dist_map;
    for (it=nodes.begin();it!=nodes.end();it++){
        n_i = *it;
        for (std::size_t i=0; i<distances.size(); i++){
            total_distance += distances(n_i,i);
            // here goes map (n_i: total_distance)
        }
        node_dist_map.insert({n_i, total_distance});
        total_distance = 0;
    }
    return node_dist_map;
}


Path::Path(std::map<std::size_t,double> &cost_map, UpperTriMat &distances, std::vector<std::pair<std::size_t,std::size_t>> &vertices) : 
    is_root_(true),
    cost_(0),
    id_(-1),
    cost_map_(&cost_map),
    distances_(&distances),
    vertices_(&vertices){}

double Path::getCost() { return cost_; }

bool Path::isRoot() { return is_root_;}

bool Path::targetReached (){
    return (getId() == -2) ? true : false;
}

std::size_t Path::getId() {return id_; }

Path * Path::getParent() {return parent_;}

std::map<std::size_t, double> * Path::getCostMap() { return  cost_map_;}

std::size_t Path::getNumNodes() { return num_nodes_;}

UpperTriMat * Path::getDistances() { return distances_;}

std::vector<std::pair<std::size_t,std::size_t>> * Path::getVertices() { return vertices_; }

double Path::calculateCost(){
    bool root_reached{false};
    double cost{0};
    Path * current_path;
    if (isRoot()){
        return getCost();
    }
    current_path = parent_;
    cost = current_path->getCost();
    if  (targetReached()) {
        return cost;
    }
    double cost_new_node = cost_map_->find(id_)->second;
    std::size_t current_id;
    while (!current_path->isRoot()){
        current_id = current_path->getId(); 
        cost_new_node -= distances_->operator()(current_id, id_);
        current_path = current_path->getParent();
    }
    cost += cost_new_node;
    return cost;
}


Path::Path(Path &parent, std::size_t node): id_(node), parent_(&parent) {
    cost_map_ = parent.getCostMap();
    distances_ = parent.getDistances();
    cost_ = calculateCost();
    vertices_ = parent.getVertices();
    num_nodes_ = parent_->getNumNodes() + 1; 
}

Path::Path(Path *parent, std::size_t node): id_(node), parent_(parent) {
    cost_map_ = parent->getCostMap();
    distances_ = parent->getDistances();
    cost_ = calculateCost();
    vertices_ = parent->getVertices();
    num_nodes_ = parent->getNumNodes() + 1; 
}

void Path::expand(std::list<Path> &paths, std::list<Path*> &unexplored_paths) {
    std::set<std::size_t> nodes_in_path;
    nodes_in_path.insert(getId());
    Path * current_path;
    current_path = this;
    while (!current_path->isRoot()){
        nodes_in_path.insert(current_path->getId());
        current_path = current_path->getParent();
    }
    bool finished{true};
    for (auto it=vertices_->begin(); it!=vertices_->end(); it++){
        if (!nodes_in_path.count(it->first) && !nodes_in_path.count(it->second)) {
            paths.push_front(Path(this, it->first));
            unexplored_paths.push_back(&paths.front());
            paths.push_back(Path(this, it->second));
            unexplored_paths.push_back(&paths.front());
            finished = false;
        }
    }
    if (finished){
        paths.push_back(Path(this, -2));
    }
}

void getStartNodes(std::list<Path> &paths, std::list<Path*> &unexplored_paths, UpperTriMat &distances, std::vector<std::pair<std::size_t,std::size_t>> &vertices, std::map<std::size_t,double> &cost_map){
    paths.clear();
    paths.push_back(Path(cost_map, distances, vertices));
    unexplored_paths.push_back(&paths.front());
}

bool lessThanPath::operator()(Path &p1, Path &p2){
    double cost1, cost2;
    cost1 = p1.getCost();
    cost2 = p2.getCost();
    if (cost1 < cost2){
        return true;
    } else if (cost1 == cost2){
        if (p1.getNumNodes() < p2.getNumNodes()){
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

bool lessThanPath::operator()(Path *p1, Path *p2){
    double cost1, cost2;
    cost1 = p1->getCost();
    cost2 = p2->getCost();
    if (cost1 < cost2){
        return true;
    } else if (cost1 == cost2){
        if (p1->getNumNodes() < p2->getNumNodes()){
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

std::vector<std::size_t> removeNodes(const std::vector<std::vector<double>> &points, double d_min){
    UpperTriMat distances;
    distances = getDistanceMatrix(points);
    std::map<std::size_t,double> nodes_dist_map;
    std::vector<std::pair<std::size_t,std::size_t>> vertices_to_remove;
    vertices_to_remove = filterVertices(distances, d_min);
    std::set<std::size_t> unique_nodes;
    unique_nodes = getUniqueNodes(vertices_to_remove);
    nodes_dist_map = getNodeTotalDistanceMap(unique_nodes, distances);
    std::list<Path> paths; 
    std::list<Path*> unexplored_paths; 
    getStartNodes(paths, unexplored_paths, distances, vertices_to_remove, nodes_dist_map);
    std::vector<std::size_t> nodes_to_remove;
    Path * current_path;
    current_path = unexplored_paths.front();
    std::cout << current_path->getId() << std::endl;
    current_path->expand(paths, unexplored_paths);
    unexplored_paths.pop_front();
    bool finished{false};
    while (!finished) {
        unexplored_paths.sort(lessThanPath());
        current_path = unexplored_paths.front();
        std::cout << current_path->getCost() << std::endl;
        if (current_path->targetReached()){
            finished = true;
        } else {
            current_path->expand(paths, unexplored_paths);
            unexplored_paths.pop_front();
        } 
    }
    return nodes_to_remove;
}

void greedyNodeRemoval(std::vector<std::vector<double>> &samples, double d_min, std::vector<bool> &accepted, int &num_accepted){
    std::set<std::size_t> removed;
    num_accepted = 0;
    double d_ij;
    int choice;
    std::size_t n = samples.size();
    for (std::size_t i=0; i<n; i++){
        for (std::size_t j=0; j<n; j++){
            if (i>j){
                if (removed.count(i) == 0 && removed.count(j) == 0) {
                    d_ij = distance(samples[i],samples[j]);
                    if (d_ij < d_min) {
                        choice = rand() % 2;
                        if (choice == 0) {
                            removed.insert(i);
                        } else {
                            removed.insert(j);
                        }
                    }
                }
            }
        }
    }
    std::size_t j{0};
    accepted.resize(samples.size());
    for (std::size_t i=0; i<n; i++){
        if (removed.count(i) == 0){
            accepted[i] = true;
            num_accepted += 1;
        }
    }
}




