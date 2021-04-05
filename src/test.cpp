#include <Eigen/Dense>
#include <vector>
#include <iostream>

using namespace Eigen;


struct ad{
    double d{0};
};

void func(void *data){
    ad* b = (ad*) data;
    std::cout << b->d <<std::endl;
    b->d = 2;
    std::cout << b->d <<std::endl;
}

int main()
{
    std::vector<double> a{1,2};
    Vector2d b;
    for (std::size_t i = 0; i < a.size(); ++i){
        b(i) = a[i];
    };
    std::cout << b << std::endl;

    ad da;
    func(&da);
    std::cout << da.d << std::endl;
}
