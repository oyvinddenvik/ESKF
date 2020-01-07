#pragma once
#include<Eigen/Dense>

using namespace Eigen;

template <typename Derived>
MatrixXd blkdiag(const MatrixBase<Derived>& a, int count)
{
    MatrixXd bdm = MatrixXd::Zero(a.rows() * count, a.cols() * count);
    for (int i = 0; i < count; ++i)
    {
        bdm.block(i * a.rows(), i * a.cols(), a.rows(), a.cols()) = a;
    }

    return bdm;
}



