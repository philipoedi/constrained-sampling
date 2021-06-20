#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "tangent.hpp"

int main(){
    using namespace Eigen;
    using namespace std;
    
    const int m1{2};
    const int m2{1};
    const int n1{3};
    const int n2{4};

    Matrix<double,m1,n1> jac1 =  Matrix<double,m1,n1>::Constant(1);
    cout << "jac1: \n" << jac1 << endl;
    Matrix<double,m2,n2> jac2 = Matrix<double,m2,n2>::Constant(1); 

    Matrix<double,n1*(n1-m1),1> theta1;
    Matrix<double,n2*(n2-m2),1> theta2;

    TangentConstraintData<n1,m1> d1;
    TangentConstraintData<n2,m2> d2;
    cout << "d1.A :\n" << d1.A << "\n d1.b : \n" << d1.b << endl;
    cout << "d2.A :\n" << d2.A << "\n d2.b : \n" << d2.b << endl;
    cout << "tangent space" << endl;
    
    cout << theta1 << endl;
    // TangentSpace
    TangentSpace<n1,m1> tp1;
    TangentSpace<n2,m2> tp2;

    tp1.createConstraintData(jac1);
    cout << "tp1: \n "<<tp1.getConstraintData().A << endl;
    tp2.createConstraintData(jac2);
    cout << "tp2: \n" << tp2.getConstraintData().A << endl;

    Matrix<double,n1,(n1-m1)> theta1m = Matrix<double,n1,(n1-m1)>::Constant(2);
    tp1.updateConstraintData(theta1m);
    cout << "Updated tp1: \n" << tp1.getConstraintData().A << endl;

    Matrix<double,n2,(n2-m2)> theta2m = Matrix<double,n2,(n2-m2)>::Constant(2);
    tp2.updateConstraintData(theta2m);
    cout << "Updated tp1: \n" << tp2.getConstraintData().A << endl;
    
    cout << d1.A*theta1 <<endl;

    cout<< "Tp1 slack: " << tp1.getSlack(tp1.getConstraintData().A,theta1, tp1.getConstraintData().b) << endl;
    
    cout<< "Tp1 gradein: \n" << tp1.getGradient(tp1.getConstraintData().A) << endl;
    cout<< "Tp2 gradein: \n" << tp2.getGradient(tp2.getConstraintData().A) << endl;

    for(int i=0; i<tp1.getConstraintData().A.size(); i++){
        cout << tp1.getConstraintData().A(i) << endl;
    }


    Map< const VectorXd> v1(tp1.getConstraintData().A.transpose().data(), tp1.getConstraintData().A.size());
    cout << v1 << endl;

    double result1[n1*(n1-m1)];
    double grad[((n1-m1)*m1+(n1-m1)*(n1-m1))*n1];
    double x[(n1-m1)*m1+(n1-m1)*(n1-m1)];
    tangentSpaceConstraints<n1,m1>(2, result1,2,x,grad,&tp1);
    
    

    Matrix<double,1,2> jac_circle;
    jac_circle << 0,2;
    
    const int n_circle{2};
    const int m_circle{1};
    TangentSpace<2,1> tpc;

    Vector2d x0{0,1};
    tpc.findTangentSpace(x0, jac_circle);
    //tp1.orthonormalBasisConstraints();

    Matrix<double,1,1> u1{1};
    cout << "map from tangent to ambient: " << endl;
    cout << tpc.toAmbient(u1) << endl;
    //

   /* 
    // createTangentConstraintData
    tp1.createConstraintData(jac1);
    tp2.createConstraintData(jac2);
        
    TangentConstraintData<> d1;
    TangentConstriantData<> d2;

    d1 = tp1.getConstraintData();
    d2 = tp2.getConstraintData();
    
    cout <<"d1.A: " << d1.A << "d1.b: " << d1.b << endl;
    cout <<"d2.A: " << d2.A << "d2.b: " << d2.b << endl;
   */
    
    // updateTangentConstraintData





};
