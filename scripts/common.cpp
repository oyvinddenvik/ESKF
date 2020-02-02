#include"common.h"







Matrix3d crossProductMatrix(const Vector3d& n)
{
	Matrix3d skewMatrix = Matrix3d::Zero();

	skewMatrix << 0, -n(2), n(1),
				  n(2), 0, -n(0),
				 -n(1), n(0), 0;

	return skewMatrix;
}

Vector4d quaternionHamiltonProduct(VectorXd quatLeft, VectorXd quatRight)
{
	//int realPartIndex = 0;

	Vector4d quatProduct = Vector4d::Zero();
	Vector3d imaginaryPartOfLeftQuaternion = Vector3d::Zero();
	Vector3d imaginaryPartOfRightQuaternion = Vector3d::Zero();
	

	if (quatLeft.size() == 3) // Assume pure quaternion
	{
		quatLeft << 0,
			quatLeft(0),
			quatLeft(1),
			quatLeft(2);
	}
	if (quatRight.size() == 3)  // Assume pure quaternion
	{
		quatRight << 0,
			quatRight(0),
			quatRight(1),
			quatRight(2);
	}
	
	imaginaryPartOfLeftQuaternion = quatLeft.block<3, 1>(1, 0);
	imaginaryPartOfRightQuaternion = quatRight.block<3, 1>(1, 0);

	
	quatProduct << quatLeft(0) * quatRight(0) - imaginaryPartOfLeftQuaternion.transpose() * imaginaryPartOfRightQuaternion,
				   quatRight(0)* imaginaryPartOfLeftQuaternion + quatLeft(0) * imaginaryPartOfRightQuaternion + crossProductMatrix(imaginaryPartOfLeftQuaternion) * imaginaryPartOfRightQuaternion;
					
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

Matrix3d quaternion2Rotationmatrix(const Vector4d& quaternion)
{
	
	Matrix3d RotationMatrix = Matrix3d::Zero();
	Matrix3d S = Matrix3d::Zero();
	Vector3d imaginaryPart = Vector3d::Zero();
	double realPart{ 0 };


	realPart = quaternion(0);
	imaginaryPart = quaternion.segment(1, quaternion.size() - 1);
	

	S = crossProductMatrix(imaginaryPart);
	//std::cout << S << std::endl;
	RotationMatrix = MatrixXd::Identity(3, 3) + (2 * realPart * S) + (2 * S * S);
	//std::cout << RotationMatrix<< std::endl;

	return RotationMatrix;

}

MatrixXd jacobianFdOfDVL(const VectorXd& fun,const VectorXd& x, const double& step, const Vector3d& velWorld )
{
	Vector4d xPerturbed = Vector4d::Zero();
	Vector3d fPerturbed = Vector3d::Zero();
	int numberOfRowsOfQuaternion{ 0 };
	int numberOfRowsOfF{ 0 };

	numberOfRowsOfQuaternion = x.rows();
	numberOfRowsOfF = fun.rows();

	MatrixXd jacobianMatrix(numberOfRowsOfF, numberOfRowsOfQuaternion);

	
	for (int i = 0; i < numberOfRowsOfQuaternion; ++i)
	{
		xPerturbed = x;
		xPerturbed(i) = xPerturbed(i) + step;
		fPerturbed = quaternion2Rotationmatrix(xPerturbed) * velWorld;
		jacobianMatrix.col(i) = (fPerturbed - fun) / step; // Check if error
	}

	
	return jacobianMatrix;

}



