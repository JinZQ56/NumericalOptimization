#include <iostream>
#include <random>
#include <cmath>
#include <Eigen/Eigen>

void rand_permutation(const int n, int *p)
{
    typedef std::uniform_int_distribution<int> rand_int;
    typedef rand_int::param_type rand_range;
    static std::mt19937_64 gen;
    static rand_int rdi(0, 1);
    int j, k;
    for (int i = 0; i < n; ++i)
    {
        p[i] = i;
    }
    for (int i = 0; i < n; ++i)
    {
        rdi.param(rand_range(0, n - i - 1));
        j = rdi(gen) + i;
        k = p[j];
        p[j] = p[i];
        p[i] = k;
    }
}

int main()
{
    int n=7;
    Eigen::VectorXi perm(n - 1);
    Eigen::VectorXi next(n);
    Eigen::VectorXi prev(n + 1);
    if (n > 1)
    {
        rand_permutation(n - 1, perm.data());
        prev(0) = 0;
        next(0) = perm(0) + 1;
        prev(perm(0) + 1) = 0;
        for (int i = 0; i < n - 2; ++i)
        {
            next(perm(i) + 1) = perm(i + 1) + 1;
            prev(perm(i + 1) + 1) = perm(i) + 1;
        }
        next(perm(n - 2) + 1) = n;
    }
    else
    {
        prev(0) = 0;
        next(0) = 1;
        next(1) = 1;
    }

    std::cout << "prem" << std::endl;
    std::cout << perm.transpose() << std::endl;

    std::cout << "next" << std::endl;
    std::cout << next.transpose() << std::endl;

    std::cout << "prev" << std::endl;
    std::cout << prev.transpose() << std::endl;
    

    return 0;

}