#include <iostream>
#include <functional>

using namespace std;

class foo
{

    public:
        foo(){};
        double add(double y) { return y + x;};

    private:
        double x{2};
};

int main()
{
    foo k;
    cout << k.add(3) << endl;
    std::function<double(double)> f = std::bind(k::add, this);
    cout << f(2) << endl;
    return 0;
};
