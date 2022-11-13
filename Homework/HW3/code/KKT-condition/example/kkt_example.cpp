#include <iostream>
#include "kkt.hpp"

using namespace std;
using namespace kkt;
using namespace Eigen;

int main(int argc, char **argv)
{
    const int d=3;
    const int n=2;
    Matrix<double, d, d> Q;
    Matrix<double, d, 1> c;
    Matrix<double, n, d> A;
    Matrix<double, n, 1> b;
    Matrix<double, d, 1> x;

    Q << 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0;
    c << 1.2, 2.5, 3.8;

    A << 1.0, 2.0, 3.0,
         2.0, 3.0, 1.0;
    b << 4.0, 5.0;

    double min = solve<d, n>(Q, c, A, b, x);
    
    cout << "optimal sol: " << x.transpose() << endl;
    cout << "optimal obj: " << min << endl;
    cout << "cons precision: " << (A * x - b).maxCoeff() << endl;

    return 0;
}