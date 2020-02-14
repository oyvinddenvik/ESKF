#pragma once
#include<Eigen/Dense>
#include <iostream>
#include <math.h>

constexpr long double PI{3.1415926535};

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


MatrixXd blk3x3Diag(const Matrix3d& matrixA, const Matrix3d& matrixB, const Matrix3d& matrixC, const Matrix3d& matrixD);

Matrix3d crossProductMatrix(const Vector3d& n);
Vector4d quaternionHamiltonProduct(VectorXd quatLeft, VectorXd quatRight);
Matrix3d quaternion2Rotationmatrix(const Vector4d& quaternion);
MatrixXd jacobianFdOfDVL(const VectorXd& fun, const VectorXd& x, const double& step, const Vector3d& velWorld);
Matrix3d eulerToRotationMatrix(const Vector3d& eulerAngles);
