#ifndef KKT_HPP
#define KKT_HPP

#include <Eigen/Eigen>

namespace kkt
{
    template<int d, int n>
    inline double solve(const Eigen::Matrix<double, d, d> &Q,
                        const Eigen::Matrix<double, d, 1> &c,
                        const Eigen::Matrix<double, n, d> &A,
                        const Eigen::Matrix<double, n, 1> &b,
                        Eigen::Matrix<double, d, 1> &x)
    {
        Eigen::Matrix<double, d+n, d+n> QA;
        Eigen::Matrix<double, d+n, 1> cb;
        Eigen::Matrix<double, d+n, 1> xv;
        double minimum;

        QA.setZero();
        QA.block(0, 0, d, d) = Q;
        QA.block(0, d, d, n) = A.transpose();
        QA.block(d, 0, n, d) = A;

        cb.setZero();
        cb.block(0, 0, d, 1) = -c;
        cb.block(d, 0, n, 1) = b;

        xv = QA.inverse() * cb;
        x = xv.block(0, 0, d, 1);

        minimum = x.transpose() * Q * x;
        minimum = 0.5 * minimum + c.transpose() * x;

        return minimum;  
    }
                 
}

#endif