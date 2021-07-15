#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <list>
#include <cassert>
#include "filter.hpp"

using namespace Eigen;
using namespace std;



int main(){
    size_t n{6};
    UpperTriMat utm;
    cout << "init utm2" << endl;
    UpperTriMat utm2(n);
    cout << "resize" << endl;
    utm.resize(n);
    cout << "Size utm: " << utm.size() << endl;
    cout << "Size utm2: " << utm.size() << endl;
    
    utm.addData(0,1,2);
    cout << "Added Data 2 at 0,1" << endl;
    cout << "1,0: " << utm(1,0) << endl;
    cout << "0,1: " << utm(0,1) << endl;

    vector<vector<double>> data;
    data.push_back(vector<double>{1,1});
    data.push_back(vector<double>{2,2});
    data.push_back(vector<double>{3,5});
    data.push_back(vector<double>{4,4});
    data.push_back(vector<double>{6,0});
    data.push_back(vector<double>{5,1});
    cout << "data[0][0]: " << data[0][1] << endl;
     
    cout << "distance between 1,1 and 2,2 " << distance(data[0],data[1]) << endl;

    utm2 = getDistanceMatrix(data);

    cout << "1,0: " << utm2(1,0) << endl;
    cout << "2,0: " << utm2(2,0) << endl;
    cout << "3,0: " << utm2(0,3) << endl;

    vector<pair<size_t,size_t>> verts;
    double d_min{2};
    vector<double> row;
    row = utm(0);

    verts = filterVertices(utm2, d_min);
    vector<pair<size_t,size_t>>::iterator it;
    pair<size_t,size_t> v;
    for (it=verts.begin(); it!=verts.end(); it++){
        v = *it;
        cout << v.first << " " << v.second << endl;
        cout << d_min << " > " << utm2(v.first,v.second) << endl;
        assert (utm(v.first, v.second) <= d_min);
    }

    set<size_t> uniques;
    uniques = getUniqueNodes(verts);
    set<size_t>::iterator it2;
    for (it2 = uniques.begin(); it2!=uniques.end(); it2++){
        cout << "unique: " << *it2 << endl;
    }
    
    map<std::size_t,double> distancemap;
    distancemap = getNodeTotalDistanceMap(uniques, utm2);
    map<size_t,double>::iterator it3;
    for (it3=distancemap.begin(); it3!=distancemap.end(); it3++){
        cout << it3->first << ": " << it3->second << endl;
    }
    
    cout << "total dist 0 = " << distancemap.find(0)->second << endl;

    UpperTriMat * utm_ptr = &utm2;
    cout << utm_ptr->operator()(0,1) << endl;


    Path p(distancemap,utm2,verts);
    
    cout << "Cost: " << p.calculateCost() << endl;

    Path p2(p,0);
    cout << "P2 Cost: " << p2.getCost() << endl;
    Path p3(p2,1);
    cout << "P3 Cost: " << p3.getCost() << endl;
    Path p4(p3,2);
    cout << "P4 Cost: " << p4.getCost() << endl;
    list<Path> expanded_paths;
    list<Path*> unexplored_paths_new;
    p2.expand(expanded_paths, unexplored_paths_new);    

    list<Path> paths_new;
    unexplored_paths_new.clear();
    getStartNodes(paths_new, unexplored_paths_new, utm2, verts, distancemap);
    Path current_path;
    current_path = paths_new.front();
    paths_new.pop_front();
    cout << "list after pop: "<< paths_new.empty() << endl; 
    cout << "Root node cost: " << current_path.getCost() << endl;
    current_path.expand(paths_new, unexplored_paths_new);
    cout << "size of updated unexplored paths list" << unexplored_paths_new.size() << endl;
    
    lessThanPath ltp;
    cout << "check less than: p2<p3 " << ltp(p2,p3) << endl;
    assert (ltp(p2,p3)) ; 
    paths_new.sort(lessThanPath());
    double c1{0}, c2;
    for (auto it=paths_new.begin();it!=paths_new.end(); it++ ){
        c2 = it->getCost();
        cout <<  it->getCost() << endl;
        assert (c2 > c1);
        c1 = c2;
    }
    
    removeNodes(data, 2);
    return 0;
}
