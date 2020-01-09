#include"common.h"


Matrix3d crossProductMatrix(Vector3d n)
{
	Matrix3d skewMatrix = Matrix3d::Zero();

	skewMatrix << 0, -n(2), n(1),
				  n(2), 0, -n(0),
				 -n(1), n(0), 0;

	return skewMatrix;
}

Vector4d quaternionHamiltonProduct(VectorXd quatLeft, VectorXd quatRight)
{
	Vector4d quatProduct = Vector4d::Zero();

	if (quatLeft.count() == 3) // Assume pure quaternion
	{
		quatLeft << 0,
			quatLeft(0),
			quatLeft(1),
			quatLeft(2);
	}
	if (quatRight.count() == 3)  // Assume pure quaternion
	{
		quatRight << 0,
			quatRight(0),
			quatRight(1),
			quatRight(2);

	}

	

	return quatProduct;



}

MatrixXd blk3x3Diag(const Matrix3d& matrixA, const Matrix3d& matrixB, const Matrix3d& matrixC, const Matrix3d& matrixD)
{
	int numberOfMatricies{ 4 };

	MatrixXd bdm = MatrixXd::Zero(matrixA.rows() * numberOfMatricies, matrixA.cols() * numberOfMatricies);

	for (int i = 0; i < numberOfMatricies; ++i)
	{
		if (i == 0)
		{
			bdm.block(i * matrixA.rows(), i * matrixA.cols(), matrixA.rows(), matrixA.cols()) = matrixA;
		}
		else if (i == 1)
		{
			bdm.block(i * matrixB.rows(), i * matrixB.cols(), matrixB.rows(), matrixB.cols()) = matrixB;
		}
		else if (i == 2)
		{
			bdm.block(i * matrixC.rows(), i * matrixC.cols(), matrixC.rows(), matrixC.cols()) = matrixC;
		}
		else if (i == 3)
		{
			bdm.block(i * matrixD.rows(), i * matrixD.cols(), matrixD.rows(), matrixD.cols()) = matrixD;
		}
		else
		{
			std::cout << "ERROR: Could not make a blkdiag matrix" << std::endl;
		}

		//bdm.block(i * a.rows(), i * a.cols(), a.rows(), a.cols()) = a;
	}
	

	return bdm;
}





