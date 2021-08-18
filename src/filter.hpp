#ifndef FILTER_H
#define FILTER_H

#include <iostream>
#include <Eigen/Dense>
#include <list>
#include <vector>
#include <cassert>
#include <set>
#include <map>
#include <cstdlib>


using namespace Eigen;


class UpperTriMat {

    public:

        UpperTriMat(){;};
        UpperTriMat(std::size_t n);
        void resize(std::size_t n);
        void addData(std::size_t i, std::size_t j, double x);
        double operator()(std::size_t i, std::size_t j);  
        std::vector<double> operator()(std::size_t i);
        std::size_t size();

    private:

        std::vector<std::vector<double>> data_;
};

class Path {

    public:
        
        Path(){};
        Path(std::map<std::size_t,double> &cost_map, UpperTriMat &d, std::vector<std::pair<std::size_t,std::size_t>> &vertices);
        Path(Path &p, std::size_t node); 
        Path(Path *p, std::size_t node); 
        double getCost();
        double calculateCost();
        bool isRoot();
        std::size_t getId();
        Path * getParent();
        std::map<std::size_t,double>* getCostMap();
        std::size_t getNumNodes();
        UpperTriMat * getDistances();
        std::vector<std::pair<std::size_t,std::size_t>> * getVertices();
        void expand(std::list<Path> &paths, std::list<Path*> &unexplored_paths);
        bool targetReached();

    private:

        std::size_t id_;
        double cost_;
        Path * parent_;
        UpperTriMat * distances_;
        std::map<std::size_t, double> * cost_map_;
        std::size_t num_nodes_;
        std::vector<std::pair<std::size_t,std::size_t>> * vertices_;
        bool is_root_{false};
};

struct lessThanPath {
    bool operator()(Path &p1,  Path &p2);
    bool operator()(Path *p1,  Path *p2);
};

double distance(const std::vector<double> &a, const std::vector<double> &b);

UpperTriMat getDistanceMatrix(const std::vector<std::vector<double>> &samples);

std::vector<std::pair<std::size_t,std::size_t>> filterVertices(UpperTriMat distances, double d_min);

std::set<std::size_t> getUniqueNodes(std::vector<std::pair<std::size_t,std::size_t>> &vertices);

std::map<std::size_t,double> getNodeTotalDistanceMap(std::set<std::size_t> &nodes, UpperTriMat &distances);
    
void getStartNodes(std::list<Path> &paths, std::list<Path*> &unexplored_paths, UpperTriMat &distances, std::vector<std::pair<std::size_t,std::size_t>> &vertices, std::map<std::size_t,double> &cost_map);

std::vector<std::size_t> removeNodes(const std::vector<std::vector<double>> &points, double d_min);

void greedyNodeRemoval(std::vector<std::vector<double>> &samples, double d_min, std::vector<bool> &accepted, int &num_accepted);

/*
class Path {

    public:
        
        Path();
        Path(std::size_t node);
        Path(Path &p, std::size_t node); 
        expand();
        getCost();
        getId();

    private:

        std::size_t id_;
        double cost_;
        Path parent_;
        std::size_t num_nodes_;

}
*/

// Method Declarations
// ----------------------------------

// others
// ----------------------------------


// rename
/*

void getStartNodes(std::list<Node> &nodes, std::map<std::size_t, double> &nodes_dist_map){
    std::set<std::size_t>>::iterator it;
    for (it=nodes_dist_map.begin(); it!= nodes_dist_map.end();it++) {
        nodes.push_back(Node(*it->first,*it->second));
}

std::vector<std::size_t> removeNodes(std::vector<std::vector<double>> &points, double d_min ){



*/

#endif
