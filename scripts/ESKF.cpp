#include "ESKF.h"
#include "common.h"



ESKF::ESKF() : pgyroBias{0},paccBias{0},SpressureZ{0},poseStates(NOMINAL_STATE_SIZE),errorStateCovariance(ERROR_STATE_SIZE,ERROR_STATE_SIZE)
{
	
	Racc.setZero();
    RaccBias.setZero();
    Rgyro.setZero();
    RgyroBias.setZero();

    Sa.setZero();
    Sg.setZero();
    Sdvl.setZero();
    Sinc.setZero();
   
	D = blk3x3Diag(Racc, Rgyro, RaccBias, RgyroBias);

	poseStates.setZero();
	errorStateCovariance.setZero();

}


ESKF::ESKF(const Matrix3d& Racc, const Matrix3d& RaccBias, const Matrix3d& Rgyro, const Matrix3d& RgyroBias, const double& pgyroBias, const double& paccBias, const Matrix3d& Sa, const Matrix3d& Sg, const Matrix3d& Sdvl, const Matrix3d& Sinc)
	:Racc{ Racc }, RaccBias{ RaccBias }, Rgyro{ Rgyro }, RgyroBias{ RgyroBias }, pgyroBias{ pgyroBias }, paccBias{ paccBias }, Sa{ Sa }, Sg{ Sg }, Sdvl{ Sdvl}, poseStates(NOMINAL_STATE_SIZE), errorStateCovariance(ERROR_STATE_SIZE,ERROR_STATE_SIZE)
{
	D = blk3x3Diag(Racc, Rgyro, RaccBias, RgyroBias);
	poseStates.setZero();
	errorStateCovariance.setZero();
	Matrix3d I_3x3 = Matrix3d::Identity();
	Vector3d initialPosition = Vector3d::Zero();
	Vector3d initialVelocity = Vector3d::Zero();
	Vector4d initialQuat = Vector4d::Zero();
	Vector3d initialAccBias = Vector3d::Zero();
	Vector3d initialGyroBias = Vector3d::Zero();
	Matrix3d initialPPos = Matrix3d::Zero();
	Matrix3d initialPVel = Matrix3d::Zero();
	Matrix3d initialPDq = Matrix3d::Zero();
	Matrix3d initialPAccBias = Matrix3d::Zero();
	Matrix3d initialPGyroBias = Matrix3d::Zero();

	initialPosition << -1.96,0.061,0;
	initialVelocity << 0,0,0;
	initialQuat << 0.9935,0,0,-0.1134;
	initialAccBias << 0,0,0;
	initialGyroBias << 0,0,0;
	initialPPos.diagonal() << 1e-9*1e-9,1e-9*1e-9,1e-9*1e-9;
	initialPVel.diagonal() << 1e-4*1e-4,1e-4*1e-4,1e-4*1e-4;
	initialPDq.diagonal() << 12e-3*12e-3,12e-3*12e-3,12e-3*12e-3;
	initialPAccBias.diagonal() << 12e-9*12e-9,12e-9*12e-9,12e-9*12e-9;
	initialPGyroBias.diagonal() << 3e-18*3e-18,3e-18*3e-18,3e-18*3e-18;
	
	errorStateCovariance.block<3,3>(0,0) = initialPPos;
	errorStateCovariance.block<3,3>(3,3) = initialPVel;
	errorStateCovariance.block<3,3>(6,6) = initialPDq;
	errorStateCovariance.block<3,3>(9,9) = initialPAccBias;
	errorStateCovariance.block<3,3>(12,12) = initialPGyroBias;

	poseStates << initialPosition,
				  initialVelocity,
				  initialQuat,
				  initialAccBias,
				  initialGyroBias;	

	//std::cout<<poseStates<<std::endl;
	//std::cout<<errorStateCovariance<<std::endl;
}

VectorXd ESKF::predictNominal(const VectorXd& xnominal, const Vector3d& accRectifiedMeasurements, const Vector3d& gyroRectifiedmeasurements, const double& Ts)
{
	
	// Initilize
	VectorXd xNextnominal(NOMINAL_STATE_SIZE); 
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
	position = xnominal.block<NOMINAL_POSITION_STATE_SIZE, 1>(NOMINAL_POSITION_STATE_OFFSET, 0);
	velocity = xnominal.block<NOMINAL_VELOCITY_STATE_SIZE, 1>(NOMINAL_VELOCITY_STATE_OFFSET, 0);
	quaternion = xnominal.block<NOMINAL_QUATERNION_STATE_SIZE, 1>(NOMINAL_QUATERNION_STATE_OFFSET, 0);
	accBias = xnominal.block<NOMINAL_ACC_BIAS_SIZE, 1>(NOMINAL_ACC_BIAS_STATE_OFFSET, 0);
	gyroBias = xnominal.block<NOMINAL_GYRO_BIAS_SIZE, 1>(NOMINAL_GYRO_BIAS_STATE_OFFSET, 0);


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

MatrixXd ESKF::Aerr(const VectorXd& xnominal, const Vector3d& accRectifiedMeasurements, const Vector3d& gyroRectifiedmeasurements)
{
	// Initilize
	MatrixXd A(ERROR_STATE_SIZE,ERROR_STATE_SIZE);
	A.setZero();
	Matrix3d rotationMatrix = Matrix3d::Zero();
	Vector4d quaternion = Vector4d::Zero();
	Matrix3d identityMatrix = Matrix3d::Identity();


	quaternion = xnominal.block<NOMINAL_QUATERNION_STATE_SIZE, 1>(NOMINAL_QUATERNION_STATE_OFFSET, 0);
	rotationMatrix = quaternion2Rotationmatrix(quaternion);


	A.block<3, 3>(0, 3) = identityMatrix;
	A.block<3, 3>(3, 6) = -1.0 * rotationMatrix * crossProductMatrix(accRectifiedMeasurements);
	A.block<3, 3>(3, 9) = -1.0 * rotationMatrix;
	A.block<3, 3>(6, 6) = -1.0 * crossProductMatrix(gyroRectifiedmeasurements);
	A.block<3, 3>(6, 12) = -identityMatrix;
	A.block<3, 3>(9, 9) = -1.0 * paccBias * identityMatrix;
	A.block<3, 3>(12, 12) = -1.0 * pgyroBias * identityMatrix;

	// Bias corrections

	A.block<3, 3>(3, 9) = A.block<3, 3>(3, 9) * Sa;
	A.block<3, 3>(6, 12) = A.block<3, 3>(6, 12) * Sg;

	return A;
} 

MatrixXd ESKF::Gerr(const VectorXd& xnominal)
{
	// Initilize
	MatrixXd Gerror(ERROR_STATE_SIZE, 12);
	Gerror.setZero();
	Matrix3d rotationMatrix = Matrix3d::Zero();
	Vector4d quaternion = Vector4d::Zero();
	Matrix3d identityMatrix = Matrix3d::Identity();


	quaternion = xnominal.block<NOMINAL_QUATERNION_STATE_SIZE, 1>(NOMINAL_QUATERNION_STATE_OFFSET, 0);
	rotationMatrix = quaternion2Rotationmatrix(quaternion);

	Gerror.block<3, 3>(3, 0) = -1.0 * rotationMatrix;
	Gerror.block<3, 3>(6, 3) = -1.0 * identityMatrix;
	Gerror.block<3, 3>(9, 6) = identityMatrix;
	Gerror.block<3, 3>(12, 9) = identityMatrix;

	return Gerror;

}

AdandGQGD ESKF::discreteErrorMatrix(const VectorXd& xnominal,const Vector3d& accRectifiedMeasurements,const Vector3d& gyroRectifiedmeasurements,const double& Ts)
{
	// Initilize
	AdandGQGD errorMatrix;
	MatrixXd vanLoan(30, 30);
	MatrixXd vanLoanExponentional(30, 30);
	MatrixXd zeros(ERROR_STATE_SIZE, ERROR_STATE_SIZE);
	MatrixXd A(ERROR_STATE_SIZE, ERROR_STATE_SIZE);
	MatrixXd G(ERROR_STATE_SIZE, 12);

	A.setZero();
	G.setZero();
	vanLoanExponentional.setZero();
	zeros.setZero();
	vanLoan.setZero();
	errorMatrix.Ad.setZero();
	errorMatrix.GQGD.setZero();

	// Caculate Van Loan
	A = Aerr(xnominal, accRectifiedMeasurements, gyroRectifiedmeasurements);
	G = Gerr(xnominal);

	vanLoan << -1.0 * A, G* D* G.transpose(),
				zeros, A.transpose();

	vanLoan = vanLoan * Ts;
	vanLoanExponentional= vanLoan.exp();
	errorMatrix.Ad = vanLoanExponentional.block<ERROR_STATE_SIZE, ERROR_STATE_SIZE>(ERROR_STATE_SIZE, ERROR_STATE_SIZE).transpose();
	errorMatrix.GQGD = vanLoanExponentional.block<ERROR_STATE_SIZE, ERROR_STATE_SIZE>(ERROR_STATE_SIZE, ERROR_STATE_SIZE).transpose() * vanLoanExponentional.block<ERROR_STATE_SIZE, ERROR_STATE_SIZE>(0, ERROR_STATE_SIZE);

	return errorMatrix;
}

MatrixXd ESKF::predictCovariance(const VectorXd& xnominal,const MatrixXd& P,const Vector3d& accRectifiedMeasurements,const Vector3d& gyroRectifiedmeasurements,const double& Ts)
{
	MatrixXd Pprediction(ERROR_STATE_SIZE, ERROR_STATE_SIZE);
	AdandGQGD errorMatrix;

	Pprediction.setZero();
	errorMatrix.Ad.setZero();
	errorMatrix.GQGD.setZero();

	errorMatrix = discreteErrorMatrix(xnominal, accRectifiedMeasurements, gyroRectifiedmeasurements, Ts);
	Pprediction = (errorMatrix.Ad * P * errorMatrix.Ad.transpose()) + errorMatrix.GQGD;

	return Pprediction;

}

void ESKF::predict(Vector3d zAccMeasurements, Vector3d zGyroMeasurements,const double& Ts)
{
	
	StatesAndErrorCovariance predictions;
	predictions.P.setZero();
	predictions.X.setZero();

	Vector3d accelerationRectified = Vector3d::Zero();
	Vector3d gyroRectified = Vector3d::Zero();

	Vector3d accBias = Vector3d::Zero();
	Vector3d gyroBias = Vector3d::Zero();

	zAccMeasurements = Sa * zAccMeasurements;
	zGyroMeasurements = Sg * zGyroMeasurements;

	accBias = Sa * poseStates.block<NOMINAL_ACC_BIAS_SIZE, 1>(NOMINAL_ACC_BIAS_STATE_OFFSET, 0);
	gyroBias = Sg * poseStates.block<NOMINAL_GYRO_BIAS_SIZE, 1>(NOMINAL_GYRO_BIAS_STATE_OFFSET, 0);

	accelerationRectified = zAccMeasurements - accBias;
	gyroRectified = zGyroMeasurements - gyroBias;


	predictions.X = predictNominal(poseStates, accelerationRectified, gyroRectified, Ts);
	predictions.P = predictCovariance(poseStates, errorStateCovariance, accelerationRectified, gyroRectified, Ts);

	poseStates = predictions.X;
	errorStateCovariance = predictions.P;
}

StatesAndErrorCovariance ESKF::inject(const VectorXd& xnominal,const VectorXd& deltaX,const MatrixXd& P)
{
	MatrixXd Ginject(ERROR_STATE_SIZE, ERROR_STATE_SIZE);
	Ginject.setIdentity();
	Matrix3d identityMatrix3x3 = Matrix3d::Identity();
	Vector3d positionInjections = Vector3d::Zero();
	Vector3d velocityInjections = Vector3d::Zero();
	Vector4d quaternionInjections = Vector4d::Zero();
	Vector3d accelerationBiasInjections = Vector3d::Zero();
	Vector3d gyroBiasInjections = Vector3d::Zero();
	Vector4d quatRight = Vector4d::Zero();
	StatesAndErrorCovariance injections;
	injections.P.setZero();
	injections.X.setZero();

	quatRight << 1,
				deltaX.block<3, 1>(6, 0) / 2.0;

	
	positionInjections = xnominal.block<NOMINAL_POSITION_STATE_SIZE, 1>(NOMINAL_POSITION_STATE_OFFSET, 0) + deltaX.block<3, 1>(0, 0);
	velocityInjections = xnominal.block<NOMINAL_VELOCITY_STATE_SIZE, 1>(NOMINAL_VELOCITY_STATE_OFFSET, 0) + deltaX.block<3, 1>(3, 0);
	quaternionInjections = quaternionHamiltonProduct(xnominal.block<NOMINAL_QUATERNION_STATE_SIZE, 1>(NOMINAL_QUATERNION_STATE_OFFSET, 0), quatRight);
	accelerationBiasInjections = xnominal.block<3, 1>(10, 0) + deltaX.block<3, 1>(9, 0);
	gyroBiasInjections = xnominal.block<3, 1>(13, 0) + deltaX.block<3, 1>(12, 0);


	// Normalize quaternion
	quaternionInjections = quaternionInjections / sqrt(quaternionInjections.array().square().sum());

	
	injections.X << positionInjections,
						  velocityInjections,
						  quaternionInjections,
						  accelerationBiasInjections,
						  gyroBiasInjections;
	
	
	Ginject.block<3, 3>(6, 6) = identityMatrix3x3 - crossProductMatrix(0.5 * deltaX.block<3,1>(6, 0));
	injections.P = Ginject * P * Ginject.transpose();


	return injections;
}

InnovationPressureStates ESKF::innovationPressureZ(const VectorXd& xnominal,const MatrixXd& P,const double& zPressureZpos,const MatrixXd& RpressureZ)
{
	// Initilize
	InnovationPressureStates pressureStates;
	Matrix<double, NOMINAL_STATE_SIZE, ERROR_STATE_SIZE> X_deltaX;
	Matrix<double, 1, NOMINAL_STATE_SIZE> Hx;
	Matrix<double, 1, 2> zero1x2Matrix;
	Matrix<double, 1, 13> zero1x13Matrix;
	MatrixXd Q_deltaT(4, 3);
	MatrixXd identity6x6Matrix(6, 6);
	double zValue{ 1 };
	double eta{ 0 };
	double eps_1{ 0 };
	double eps_2{ 0 };
	double eps_3{ 0 };
	pressureStates.pressureH.setZero();
	pressureStates.pressureInnovation = 0;
	pressureStates.pressureInnovationCovariance.setZero();
	Hx.setZero();
	zero1x2Matrix.setZero();
	zero1x13Matrix.setZero();
	X_deltaX.setZero();
	identity6x6Matrix.setIdentity();
	Q_deltaT.setZero();


	eta = xnominal(6);
	eps_1 = xnominal(7);
	eps_2 = xnominal(8);
	eps_3 = xnominal(9);

	// Measurement Matrix
	Hx << zero1x2Matrix,
		  zValue,
		  zero1x13Matrix;

	X_deltaX.block<6, 6>(0, 0) = identity6x6Matrix;
	X_deltaX.block<6, 6>(10, 9) = identity6x6Matrix;
	Q_deltaT << -1.0 * eps_1, -1.0 * eps_2, -1.0 * eps_3,
				eta, -1.0 * eps_3, eps_2,
				eps_3, eta, -1.0 * eps_1,
				-1.0 * eps_2, eps_1, eta;
	X_deltaX.block<4, 3>(6, 6) = Q_deltaT;

	pressureStates.pressureH = Hx * X_deltaX;
	pressureStates.pressureInnovation = zPressureZpos - xnominal(StateMemberZ);
	pressureStates.pressureInnovationCovariance = (pressureStates.pressureH * P * pressureStates.pressureH.transpose()) + RpressureZ;


	return pressureStates;

}

void ESKF::updatePressureZ(const double& zPressureZpos,const MatrixXd& RpressureZ)
{
	StatesAndErrorCovariance injections;
	InnovationPressureStates pressureStates;
	MatrixXd identity15x15(ERROR_STATE_SIZE,ERROR_STATE_SIZE);
	MatrixXd kalmanGain(ERROR_STATE_SIZE, 1);
	MatrixXd deltaX(ERROR_STATE_SIZE, 1);
	MatrixXd pUpdate(ERROR_STATE_SIZE, ERROR_STATE_SIZE);

	identity15x15.setIdentity();
	pressureStates.pressureH.setZero();
	pressureStates.pressureInnovationCovariance.setZero();
	pressureStates.pressureInnovation = 0;
	injections.X.setZero();
	injections.P.setZero();
	kalmanGain.setZero();
	deltaX.setZero();
	pUpdate.setZero();

	pressureStates = innovationPressureZ(poseStates,errorStateCovariance, zPressureZpos, RpressureZ);
	

	// ESKF Update step
	kalmanGain = errorStateCovariance * pressureStates.pressureH.transpose() * pressureStates.pressureInnovationCovariance.inverse();
	deltaX = kalmanGain * pressureStates.pressureInnovation;
	pUpdate = (identity15x15 - (kalmanGain * pressureStates.pressureH)) * errorStateCovariance;
	injections = inject(poseStates, deltaX, pUpdate);

	poseStates = injections.X;
	errorStateCovariance = injections.P;
}

InnovationDVLStates ESKF::innovationDVL(const VectorXd& xnominal,const MatrixXd& P,const Vector3d& zDVLvel,const Matrix3d& RDVL)
{
	InnovationDVLStates dvlStates;
	double eta{ 0 };
	double eps_1{ 0 };
	double eps_2{ 0 };
	double eps_3{ 0 };
	double dx{ 0.000000001};
	Matrix3d zero3x3Matrix = Matrix3d::Zero();
	MatrixXd zero3x6Matrix(3, 6);
	Matrix<double, NOMINAL_STATE_SIZE, ERROR_STATE_SIZE> X_deltaX;
	MatrixXd Q_deltaT(4, 3);
	MatrixXd identity6x6Matrix(6, 6);

	Matrix3d R_body_to_world = Matrix3d::Zero();
	Matrix3d R_world_to_body = Matrix3d::Zero();
	Matrix3d Hv = Matrix3d::Zero();
	MatrixXd Hq(3, 4);
	Matrix<double, 3, NOMINAL_STATE_SIZE> Hx;
	Vector3d f = Vector3d::Zero();
	MatrixXd jacobianMatrix(3, 4);

	Vector4d nominalQuaternion = Vector4d::Zero();
	Vector4d q_world_to_body = Vector4d::Zero();
	Vector3d vel_world = Vector3d::Zero();
	Vector4d x = Vector4d::Zero();

	dvlStates.DVLH.setZero();
	dvlStates.DVLInnovation.setZero();
	dvlStates.DVLInnovationCovariance.setZero();
	jacobianMatrix.setZero();
	Hq.setZero();
	Hx.setZero();
	zero3x6Matrix.setZero();
	identity6x6Matrix.setIdentity();
	X_deltaX.setZero();
	Q_deltaT.setZero();

	eta = xnominal(6);
	eps_1 = xnominal(7);
	eps_2 = xnominal(8);
	eps_3 = xnominal(9);

	vel_world = xnominal.block<3, 1>(3, 0);
	nominalQuaternion = xnominal.block<4, 1>(6, 0);

	R_body_to_world = quaternion2Rotationmatrix(nominalQuaternion);
	R_world_to_body = R_body_to_world.transpose(); // Check this for error
	q_world_to_body << eta,
						-1.0 * eps_1,
						-1.0 * eps_2,
						-1.0 * eps_3;
	
	Hv = R_world_to_body;

	x = q_world_to_body;
	f = quaternion2Rotationmatrix(x) * vel_world;
	jacobianMatrix = jacobianFdOfDVL(f, x, dx,vel_world);

	Hq = jacobianMatrix;

	//std::cout << Hq << std::endl;

	Hx << zero3x3Matrix,
			Hv,
			Hq,
			zero3x6Matrix;

	//std::cout << Hx << std::endl;

	X_deltaX.block<6, 6>(0, 0) = identity6x6Matrix;
	X_deltaX.block<6, 6>(10, 9) = identity6x6Matrix;
	Q_deltaT << -1.0 * eps_1, -1.0 * eps_2, -1.0 * eps_3,
		eta, -1.0 * eps_3, eps_2,
		eps_3, eta, -1.0 * eps_1,
		-1.0 * eps_2, eps_1, eta;
	Q_deltaT = Q_deltaT * 0.5;
	X_deltaX.block<4, 3>(6, 6) = Q_deltaT;

	//std::cout << X_deltaX << std::endl;
		  
	dvlStates.DVLH = Hx * X_deltaX;
	//std::cout << dvlStates.DVLH << std::endl;
	dvlStates.DVLInnovation = zDVLvel - R_world_to_body * vel_world;
	dvlStates.DVLInnovationCovariance = (dvlStates.DVLH * P * dvlStates.DVLH.transpose()) + RDVL;

	return dvlStates;
}

void ESKF::updateDVL(const Vector3d& zDVLvel,const Matrix3d& RDVL)
{
	StatesAndErrorCovariance injections;
	InnovationDVLStates DVLstates;
	MatrixXd identity15x15(ERROR_STATE_SIZE,ERROR_STATE_SIZE);
	MatrixXd kalmanGain(ERROR_STATE_SIZE, 3);
	MatrixXd deltaX(ERROR_STATE_SIZE, 1);

	MatrixXd pUpdate(ERROR_STATE_SIZE, ERROR_STATE_SIZE);


	injections.P.setZero();
	injections.X.setZero();
	DVLstates.DVLH.setZero();
	DVLstates.DVLInnovation.setZero();
	DVLstates.DVLInnovationCovariance.setZero();
	identity15x15.setIdentity();

	DVLstates = innovationDVL(poseStates, errorStateCovariance, zDVLvel, RDVL);


	kalmanGain = errorStateCovariance * DVLstates.DVLH.transpose() * DVLstates.DVLInnovationCovariance.inverse();
	deltaX = kalmanGain * DVLstates.DVLInnovation;
	pUpdate = (identity15x15 - (kalmanGain * DVLstates.DVLH)) * errorStateCovariance;
	injections = inject(poseStates, deltaX, pUpdate);

	poseStates = injections.X;
	errorStateCovariance = injections.P;
}




