#include "ESKF.h"
#include "common.h"



ESKF::ESKF(Matrix3d Racc, Matrix3d RaccBias, Matrix3d Rgyro, Matrix3d RgyroBias, double pgyroBias, double paccBias, Matrix3d Sa, Matrix3d Sg, Matrix3d Sdvl, Matrix3d Sinc, double SpressureZ)
	:Racc{ Racc }, RaccBias{ RaccBias }, Rgyro{ Rgyro }, RgyroBias{ RgyroBias }, pgyroBias{ pgyroBias }, paccBias{ paccBias }, Sa{ Sa }, Sg{ Sg }, Sdvl{ Sdvl }, Sinc{ Sinc }, SpressureZ{SpressureZ}
{
	D = blk3x3Diag(Racc, Rgyro, RaccBias, RgyroBias);

}


VectorXd ESKF::predictNominal(VectorXd xnominal, Vector3d accRectifiedMeasurements, Vector3d gyroRectifiedmeasurements, double Ts)
{
	
	VectorXd xNextnominal(16); //(16x1)
	xNextnominal.setZero();

	Vector3d dTheta = Vector3d::Zero();
	double dAngle{0};

	Matrix3d rotationMatrix = Matrix3d::Zero();
	Vector3d position = Vector3d::Zero();
	Vector3d velocity = Vector3d::Zero();
	Vector4d quaternion = Vector4d::Zero();
	Vector3d accBias = Vector3d::Zero();
	Vector3d accGyro = Vector3d::Zero();

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
	accGyro = xnominal.block<3, 1>(13, 0);

		//std::cout << quatLeft.block<3, 1>(3, 0) << std::endl;

	// Get rotation matrix from quaternion
	rotationMatrix = quaternion2Rotationmatrix(quaternion);

	// Predictions
	predictPosition = position + (Ts * velocity) + (0.5 * Ts * Ts * ((rotationMatrix * accRectifiedMeasurements) + gravity));
	predictVelocity = velocity + Ts * ((rotationMatrix * accRectifiedMeasurements) + gravity);
	dTheta = Ts * gyroRectifiedmeasurements;
	dAngle = sqrt(dTheta.array().square().sum());



	return xNextnominal;

}


