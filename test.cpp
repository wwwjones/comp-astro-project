#include <iostream>
#include <vector>
#include <Eigen/Core>

class Val {
public:
    Val() : numbers_(5, nullptr) {}
    void set_number(int i, int* val) {numbers_[i] = val;}
    int* get_number(int i) {return numbers_[i];}
private:
    std::vector<int*> numbers_;
};

int main() {

    Val val;
    int a = 5;
    int* b = &a;
    for (int i = 0; i < 5; ++i)
        std::cout << val.get_number(i) << std::endl;

    for (int i = 0; i < 5; ++i)
        val.set_number(i, b);

    for (int i = 0; i < 5; ++i)
        std::cout << val.get_number(i) << std::endl;


    Eigen::Vector3d vec = Eigen::Vector3d::Zero();

    vec.x() = 5;

    std::cout << vec << std::endl;

    double x = 7.;
    unsigned int y = 2;

    std::cout << x/y << std::endl;

    return 0;
}
