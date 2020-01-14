#include "ESKF.h"
#include "common.h"


ESKF::ESKF(Matrix3d Racc, Matrix3d RaccBias, Matrix3d Rgyro, Matrix3d RgyroBias, double pgyroBias, double paccBias, Matrix3d Sa, Matrix3d Sg, Matrix3d Sdvl, Matrix3d Sinc, double SpressureZ)
	:Racc{ Racc }, RaccBias{ RaccBias }, Rgyro{ Rgyro }, RgyroBias{ RgyroBias }, pgyroBias{ pgyroBias }, paccBias{ paccBias }, Sa{ Sa }, Sg{ Sg }, Sdvl{ Sdvl }, Sinc{ Sinc }, SpressureZ{SpressureZ}
{
	D = blk3x3Diag(Racc, Rgyro, RaccBias, RgyroBias);
}


VectorXd ESKF::predictNominal(VectorXd xnominal, Vector3d accRectifiedMeasurements, Vector3d gyroRectifiedmeasurements, double Ts)
{
	
	// Initilize
	VectorXd xNextnominal(16); //(16x1)
	xNextnominal.setZero();

	Vector3d dTheta = Vector3d::Zero();
	double dAngle{0};
	Vector4d dq = Vector4d::Zero();

	Matrix3d rotationMatrix = Matrix3d::Zero();
	Vector3d position = Vector3d::Zero();
	Vector3d velocity = Vector3d::Zero();
	Vector4d quaternion = Vector4d::Zero();
	Vector3d accBias = Vector3d::Zero();
	Vector3d gyroBias = Vector3d::Zero();

	Vector3d predictPosition = Vector3d::Zero();
	Vector3d predictVelocity = Vector3d::Zero();
	Vector4d predictQauternion = Vector4d::Zero();
	Vector3d predictAccBias = Vector3d::Zero();
	Vector3d predictGyroBias = Vector3d::Zero();
	

	// Extract states
	position = xnominal.block<3, 1>(0, 0);
	velocity = xnominal.block<3, 1>(3, 0);
	quaternion = xnominal.block<4, 1>(6, 0);
	accBias = xnominal.block<3, 1>(10, 0);
	gyroBias = xnominal.block<3, 1>(13, 0);


	// Get rotation matrix from quaternion
	rotationMatrix = quaternion2Rotationmatrix(quaternion);

	// Predictions
	predictPosition = position + (Ts * velocity) + (0.5 * Ts * Ts * ((rotationMatrix * accRectifiedMeasurements) + gravity));
	predictVelocity = velocity + Ts * ((rotationMatrix * accRectifiedMeasurements) + gravity);
	dTheta = Ts * gyroRectifiedmeasurements;
	dAngle = sqrt(dTheta.array().square().sum());
	dq << cos(dAngle / 2.0),
		sin(dAngle / 2.0)* (dTheta / dAngle);
	predictQauternion = quaternionHamiltonProduct(quaternion, dq);
	predictAccBias = exp(-1.0*paccBias * Ts) * accBias;
	predictGyroBias = exp(-1.0*pgyroBias * Ts) * gyroBias;

	// Normalize quaternion
	predictQauternion = predictQauternion / sqrt(predictQauternion.array().square().sum());

	// Concatenate
	xNextnominal << predictPosition,
					predictVelocity,
					predictQauternion,
					predictAccBias,
					predictGyroBias;

	return xNextnominal;

}

MatrixXd ESKF::Aerr(VectorXd xnominal, Vector3d accRectifiedMeasurements, Vector3d gyroRectifiedmeasurements)
{
	// Initilize
	MatrixXd A(16,16);
	A.setZero();
	Matrix3d rotationMatrix = Matrix3d::Zero();
	Vector4d quaternion = Vector4d::Zero();
	Matrix3d identityMatrix = Matrix3d::Identity();


	quaternion = xnominal.block<4, 1>(6, 0);
	rotationMatrix = quaternion2Rotationmatrix(quaternion);


	A.block<3, 3>(0, 3) = identityMatrix;
	A.block<3, 3>(3, 6) = -1.0 * rotationMatrix * crossProductMatrix(accRectifiedMeasurements);
	A.block<3, 3>(3, 9) = -1.0 * rotationMatrix;
	A.block<3, 3>(6, 6) = -1.0 * crossProductMatrix(gyroRectifiedmeasurements);
	A.block<3, 3>(6, 12) = -1.0 * identityMatrix;
	A.block<3, 3>(9, 9) = -1.0 * paccBias * identityMatrix;
	A.block<3, 3>(12, 12) = -1.0 * pgyroBias * identityMatrix;

	// Bias corrections

	A.block<3, 3>(3, 9) = A.block<3, 3>(3, 9) * Sa;
	A.block<3, 3>(6, 12) = A.block<3, 3>(6, 12) * Sg;

	return A;
}

MatrixXd ESKF::Gerr(VectorXd xnominal)
{
	// Initilize
	MatrixXd Gerror(15, 12);
	Gerror.setZero();
	Matrix3d rotationMatrix = Matrix3d::Zero();
	Vector4d quaternion = Vector4d::Zero();
	Matrix3d identityMatrix = Matrix3d::Identity();


	quaternion = xnominal.block<4, 1>(6, 0);
	rotationMatrix = quaternion2Rotationmatrix(quaternion);

	Gerror.block<3, 3>(3, 0) = -1.0 * rotationMatrix;
	Gerror.block<3, 3>(6, 3) = -1.0 * identityMatrix;
	Gerror.block<3, 3>(9, 6) = identityMatrix;
	Gerror.block<3, 3>(12, 9) = identityMatrix;

	return Gerror;

}

AdandGQGD ESKF::discreteErrorMatrix(VectorXd xnominal, Vector3d accRectifiedMeasurements, Vector3d gyroRectifiedmeasurements, double Ts)
{
	AdandGQGD errorMatrix;
	MatrixXd vanLoan(30, 30);
	MatrixXd vanLoanExponentional(30, 30);
	MatrixXd zeros(15, 15);
	MatrixXd A(16, 16);
	MatrixXd G(15, 12);

	A.setZero();
	G.setZero();
	vanLoanExponentional.setZero();
	zeros.setZero();
	vanLoan.setZero();
	errorMatrix.Ad.setZero();
	errorMatrix.GQGD.setZero();

	A = Aerr(xnominal, accRectifiedMeasurements, gyroRectifiedmeasurements);
	G = Gerr(xnominal);

	vanLoan << -1.0 * A, G* D* G.transpose(),
				zeros, A.transpose();

	vanLoanExponentional= vanLoan.exp();
	

	errorMatrix.Ad = vanLoanExponentional.block<15, 15>(15, 15).transpose();
	errorMatrix.GQGD = vanLoanExponentional.block<15, 15>(15, 15).transpose() * vanLoanExponentional.block<15, 15>(0, 15);

	return errorMatrix;



}

MatrixXd ESKF::predictCovariance(VectorXd xnominal, MatrixXd P, Vector3d accRectifiedMeasurements, Vector3d gyroRectifiedmeasurements, double Ts)
{
	MatrixXd Pprediction(15, 15);
	AdandGQGD errorMatrix;

	Pprediction.setZero();
	errorMatrix.Ad.setZero();
	errorMatrix.GQGD.setZero();

	errorMatrix = discreteErrorMatrix(xnominal, accRectifiedMeasurements, gyroRectifiedmeasurements, Ts);
	Pprediction = (errorMatrix.Ad * P * errorMatrix.Ad.transpose()) + errorMatrix.GQGD;

	return Pprediction;

}

StatePredictions ESKF::predict(VectorXd xnominal, MatrixXd P, VectorXd zAccMeasurements, VectorXd zGyroMeasurements, double Ts)
{
	
	StatePredictions predictions;
	predictions.pPrediction.setZero();
	predictions.xNominalPrediction.setZero();

	Vector3d accelerationRectified = Vector3d::Zero();
	Vector3d gyroRectified = Vector3d::Zero();

	Vector3d accBias = Vector3d::Zero();
	Vector3d gyroBias = Vector3d::Zero();

	zAccMeasurements = Sa * zAccMeasurements;
	zGyroMeasurements = Sg * zGyroMeasurements;

	accBias = Sa * xnominal.block<3, 1>(10, 0);
	gyroBias = Sg * xnominal.block<3, 1>(13, 0);

	accelerationRectified = zAccMeasurements - accBias;
	gyroRectified = zGyroMeasurements - gyroBias;

	predictions.xNominalPrediction = predictNominal(xnominal, accelerationRectified, gyroRectified, Ts);
	predictions.pPrediction = predictCovariance(xnominal, P, accelerationRectified, gyroRectified, Ts);
	
	return predictions;
}

InjectionStates ESKF::inject(VectorXd xnominal, VectorXd deltaX, MatrixXd P)
{
	MatrixXd Ginject(15, 15);
	Ginject.setIdentity();
	Matrix3d identityMatrix3x3 = Matrix3d::Identity();
	Vector3d positionInjections = Vector3d::Zero();
	Vector3d velocityInjections = Vector3d::Zero();
	Vector4d quaternionInjections = Vector4d::Zero();
	Vector3d accelerationBiasInjections = Vector3d::Zero();
	Vector3d gyroBiasInjections = Vector3d::Zero();
	Vector4d quatRight = Vector4d::Zero();
	InjectionStates injections;
	injections.pInject.setZero();
	injections.xInject.setZero();

	quatRight << 1,
				deltaX.block<3, 1>(6, 0) / 2.0;

	
	positionInjections = xnominal.block<3, 1>(0, 0) + deltaX.block<3, 1>(0, 0);
	velocityInjections = xnominal.block<3, 1>(3, 0) + deltaX.block<3, 1>(3, 0);
	quaternionInjections = quaternionHamiltonProduct(xnominal.block<4, 1>(6, 0), quatRight);
	accelerationBiasInjections = xnominal.block<3, 1>(10, 0) + deltaX.block<3, 1>(9, 0);
	gyroBiasInjections = xnominal.block<3, 1>(13, 0) + deltaX.block<3, 1>(12, 0);


	// Normalize quaternion
	quaternionInjections = quaternionInjections / sqrt(quaternionInjections.array().square().sum());

	
	injections.xInject << positionInjections,
						velocityInjections,
						quaternionInjections,
						accelerationBiasInjections,
						gyroBiasInjections;
	
	

	Ginject.block<3, 3>(6, 6) = identityMatrix3x3 - crossProductMatrix(0.5 * deltaX.block<3,1>(6, 0));
	
	injections.pInject = Ginject * P * Ginject.transpose();

	return injections;
}

InnovationPressureStates ESKF::innovationPressureZ(VectorXd xnominal, MatrixXd P, double zPressureZpos, MatrixXd RpressureZ)
{
	InnovationPressureStates pressureStates;
	Matrix<double, 1, 16> Hx;
	double eta{ 0 };
	double eps_1{ 0 };
	double eps_2{ 0 };
	double eps_3{ 0 };
	pressureStates.pressureH.setZero();
	pressureStates.pressureInnovation.setZero();
	pressureStates.pressureInnovationCovariance.setZero();
	Hx.setZero();


	eta = xnominal(6);
	eps_1 = xnominal(7);
	eps_2 = xnominal(8);
	eps_3 = xnominal(9);

	return pressureStates;





}