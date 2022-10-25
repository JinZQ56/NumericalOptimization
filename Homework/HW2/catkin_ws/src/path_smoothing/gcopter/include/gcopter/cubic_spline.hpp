#ifndef CUBIC_SPLINE_HPP
#define CUBIC_SPLINE_HPP

#include "cubic_curve.hpp"

#include <Eigen/Eigen>

#include <cmath>
#include <vector>

namespace cubic_spline
{

    // The banded system class is used for solving
    // banded linear system Ax=b efficiently.
    // A is an N*N band matrix with lower band width lowerBw
    // and upper band width upperBw.
    // Banded LU factorization has O(N) time complexity.
    class BandedSystem
    {
    public:
        // The size of A, as well as the lower/upper
        // banded width p/q are needed
        inline void create(const int &n, const int &p, const int &q)
        {
            // In case of re-creating before destroying
            destroy();
            N = n;
            lowerBw = p;
            upperBw = q;
            int actualSize = N * (lowerBw + upperBw + 1);
            ptrData = new double[actualSize];
            std::fill_n(ptrData, actualSize, 0.0);
            return;
        }

        inline void destroy()
        {
            if (ptrData != nullptr)
            {
                delete[] ptrData;
                ptrData = nullptr;
            }
            return;
        }

    private:
        int N;
        int lowerBw;
        int upperBw;
        // Compulsory nullptr initialization here
        double *ptrData = nullptr;

    public:
        // Reset the matrix to zero
        inline void reset(void)
        {
            std::fill_n(ptrData, N * (lowerBw + upperBw + 1), 0.0);
            return;
        }

        // The band matrix is stored as suggested in "Matrix Computation"
        inline const double &operator()(const int &i, const int &j) const
        {
            return ptrData[(i - j + upperBw) * N + j];
        }

        inline double &operator()(const int &i, const int &j)
        {
            return ptrData[(i - j + upperBw) * N + j];
        }

        // This function conducts banded LU factorization in place
        // Note that NO PIVOT is applied on the matrix "A" for efficiency!!!
        inline void factorizeLU()
        {
            int iM, jM;
            double cVl;
            for (int k = 0; k <= N - 2; ++k)
            {
                iM = std::min(k + lowerBw, N - 1);
                cVl = operator()(k, k);
                for (int i = k + 1; i <= iM; ++i)
                {
                    if (operator()(i, k) != 0.0)
                    {
                        operator()(i, k) /= cVl;
                    }
                }
                jM = std::min(k + upperBw, N - 1);
                for (int j = k + 1; j <= jM; ++j)
                {
                    cVl = operator()(k, j);
                    if (cVl != 0.0)
                    {
                        for (int i = k + 1; i <= iM; ++i)
                        {
                            if (operator()(i, k) != 0.0)
                            {
                                operator()(i, j) -= operator()(i, k) * cVl;
                            }
                        }
                    }
                }
            }
            return;
        }

        // This function solves Ax=b, then stores x in b
        // The input b is required to be N*m, i.e.,
        // m vectors to be solved.
        template <typename EIGENMAT>
        inline void solve(EIGENMAT &b) const
        {
            int iM;
            for (int j = 0; j <= N - 1; ++j)
            {
                iM = std::min(j + lowerBw, N - 1);
                for (int i = j + 1; i <= iM; ++i)
                {
                    if (operator()(i, j) != 0.0)
                    {
                        b.row(i) -= operator()(i, j) * b.row(j);
                    }
                }
            }
            for (int j = N - 1; j >= 0; --j)
            {
                b.row(j) /= operator()(j, j);
                iM = std::max(0, j - upperBw);
                for (int i = iM; i <= j - 1; ++i)
                {
                    if (operator()(i, j) != 0.0)
                    {
                        b.row(i) -= operator()(i, j) * b.row(j);
                    }
                }
            }
            return;
        }

        // This function solves ATx=b, then stores x in b
        // The input b is required to be N*m, i.e.,
        // m vectors to be solved.
        template <typename EIGENMAT>
        inline void solveAdj(EIGENMAT &b) const
        {
            int iM;
            for (int j = 0; j <= N - 1; ++j)
            {
                b.row(j) /= operator()(j, j);
                iM = std::min(j + upperBw, N - 1);
                for (int i = j + 1; i <= iM; ++i)
                {
                    if (operator()(j, i) != 0.0)
                    {
                        b.row(i) -= operator()(j, i) * b.row(j);
                    }
                }
            }
            for (int j = N - 1; j >= 0; --j)
            {
                iM = std::max(0, j - lowerBw);
                for (int i = iM; i <= j - 1; ++i)
                {
                    if (operator()(j, i) != 0.0)
                    {
                        b.row(i) -= operator()(j, i) * b.row(j);
                    }
                }
            }
        }
    };

    class CubicSpline
    {
    public:
        CubicSpline() = default;
        ~CubicSpline() { A.destroy(); }

    private:
        int N;
        Eigen::Vector2d headP;
        Eigen::Vector2d tailP;
        Eigen::Matrix2Xd inP;
        BandedSystem A;
        Eigen::MatrixX2d b;
        Eigen::Matrix2Xd cMats;

    public:
        inline void setConditions(const Eigen::Vector2d &headPos,
                                  const Eigen::Vector2d &tailPos,
                                  const int &pieceNum)
        {
           //TODO
            headP = headPos;
            tailP = tailPos;
            N = pieceNum;
            A.create(N-1, 1, 1);
            b.resize(N-1, 2);
            // cMatVec.reserve(N);
            return;
        }

        inline void setInnerPoints(const Eigen::Ref<const Eigen::Matrix2Xd> &inPs)
        {
          //TODO
            inP = inPs;
            A.reset();
            b.setZero();

            A(0, 0) = 4.0;
            A(0, 1) = 1.0;
            b.row(0) = 3 * (inPs.col(1)-headP).transpose();

            for(int i=1; i<N-2; i++)
            {
                A(i, i-1) = 1.0;
                A(i, i) = 4.0;
                A(i, i+1) = 1.0;
                b.row(i) = 3 * (inPs.col(i+1)-inPs.col(i-1)).transpose();
            }

            A(N-2, N-3) = 1.0;
            A(N-2, N-2) = 4.0;
            b.row(N-2) = 3 * (tailP-inPs.col(N-3)).transpose();

            A.factorizeLU();
            A.solve(b);

            cMats.resize(2, 4*N);
            cMats.setZero();

            cMats.col(3) = headP;
            cMats.col(2) = Eigen::Vector2d::Zero();
            cMats.col(1) = 3*(inP.col(0)-headP) - b.row(0).transpose();
            cMats.col(0) = 2*(headP-inP.col(0)) + b.row(0).transpose();
            for(int i=1; i<N-1; i++)
            {
                cMats.col(4*i+3) = inP.col(i-1);
                cMats.col(4*i+2) = b.row(i-1).transpose();
                cMats.col(4*i+1) = 3*(inP.col(i)-inP.col(i-1)) - 
                                   2*b.row(i-1).transpose() - b.row(i).transpose();
                cMats.col(4*i+0) = 2*(inP.col(i-1)-inP.col(i)) + 
                                   1*b.row(i-1).transpose() + b.row(i).transpose();
            }
            cMats.col(4*(N-1)+3) = inP.col(N-2);
            cMats.col(4*(N-1)+2) = b.row(N-2).transpose();
            cMats.col(4*(N-1)+1) = 3*(tailP-inP.col(N-2)) - 2*b.row(N-2).transpose();
            cMats.col(4*(N-1)+0) = 2*(inP.col(N-2)-tailP) + 1*b.row(N-2).transpose();

            return;
        }

        inline void getCurve(CubicCurve &curve) const
        {
          //TODO
           curve.clear();
           curve.reserve(N);
           
           Eigen::Matrix<double, 2, 4> cMat;

           cMat.col(3) = headP;
           cMat.col(2) = headP - headP;
           cMat.col(1) = 3*(inP.col(0)-headP) - b.row(0).transpose();
           cMat.col(0) = 2*(headP-inP.col(0)) + b.row(0).transpose();
           curve.emplace_back(1.0, cMat);

           for(int i=1; i<N-1; i++)
           {
                cMat.col(3) = inP.col(i-1);
                cMat.col(2) = b.row(i-1).transpose();
                cMat.col(1) = 3*(inP.col(i)-inP.col(i-1)) - 
                              2*b.row(i-1).transpose() - b.row(i).transpose();
                cMat.col(0) = 2*(inP.col(i-1)-inP.col(i)) + 
                              b.row(i-1).transpose() + b.row(i).transpose();
                curve.emplace_back(1.0, cMat);
           }

           cMat.col(3) = inP.col(N-2);
           cMat.col(2) = b.row(N-2).transpose();
           cMat.col(1) = 3*(tailP-inP.col(N-2)) - 2*b.row(N-2).transpose();
           cMat.col(0) = 2*(inP.col(N-2)-tailP) + b.row(N-2).transpose();
           curve.emplace_back(1.0, cMat);
           
           return;
        }

        inline void getCMatVec(std::vector<Eigen::Matrix<double,2,4>> &cMatVec) const
        {
           Eigen::Matrix<double, 2, 4> cMat;

           cMat.col(3) = headP;
           cMat.col(2) = headP - headP;
           cMat.col(1) = 3*(inP.col(0)-headP) - b.row(0).transpose();
           cMat.col(0) = 2*(headP-inP.col(0)) + b.row(0).transpose();
           cMatVec[0] = cMat;

           for(int i=1; i<N-1; i++)
           {
                cMat.col(3) = inP.col(i-1);
                cMat.col(2) = b.row(i-1).transpose();
                cMat.col(1) = 3*(inP.col(i)-inP.col(i-1)) - 
                              2*b.row(i-1).transpose() - b.row(i).transpose();
                cMat.col(0) = 2*(inP.col(i-1)-inP.col(i)) + 
                              b.row(i-1).transpose() + b.row(i).transpose();
                cMatVec[i] = cMat;
           }

           cMat.col(3) = inP.col(N-2);
           cMat.col(2) = b.row(N-2).transpose();
           cMat.col(1) = 3*(tailP-inP.col(N-2)) - 2*b.row(N-2).transpose();
           cMat.col(0) = 2*(inP.col(N-2)-tailP) + b.row(N-2).transpose();
           cMatVec[N-1] = cMat;
        }

        inline void getStretchEnergy(double &energy) const
        {
           //TODO
            energy=0.0;
            std::vector<Eigen::Matrix<double,2,4>> cMatVec;
            cMatVec.reserve(N);
            getCMatVec(cMatVec);
            for(int i=0; i<N; i++)
            {   
                for(int j=0; j<2; j++)
                {
                    energy += 4*std::pow(cMatVec[i](j, 1), 2);
                    energy += 12*std::pow(cMatVec[i](j, 0), 2);
                    energy += 12*cMatVec[i](j,0)*cMatVec[i](j,1);
                } 
            }
            return;
        }

        inline const Eigen::MatrixX2d &getCoeffs(void) const
        {
            return b;
            // return cMats; 
        }

        inline void getGrad(Eigen::Ref<Eigen::Matrix2Xd> grad) const
        {
            //TODO
            std::vector<Eigen::Matrix<double,2,4>> cMatVec;
            cMatVec.reserve(N);
            getCMatVec(cMatVec);
            double dc0, dc1, dd0, dd1;
            for(int i=0; i<N-1; i++)
            {   
                for(int j=0; j<2; j++)
                {
                    dc0 = 8*cMatVec[i](j, 1) + 12*cMatVec[i](j, 0);
                    dd0 = 12*cMatVec[i](j, 1) + 24*cMatVec[i](j, 0);
                    dc1 = 8*cMatVec[i+1](j, 1) + 12*cMatVec[i+1](j, 0);
                    dd1 = 12*cMatVec[i+1](j, 1) + 24*cMatVec[i+1](j, 0);
                    grad(j, i) = 3*(dc0-dc1) + 2*(dd1-dd0);
                }
            }
            return;
        }
    };
}

#endif
