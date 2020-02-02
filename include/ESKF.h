#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;

enum NominalStateMembers //(16x1)
{
	StateMemberX = 0,
	StateMemberY = 1,
	StateMemberZ = 2,
	StateMemberVx = 3,
	StateMemberVy = 4,
	StateMemberVz = 5,
	StateMemberQw = 6,
	StateMemberQx = 7,
	StateMemberQy = 8,
	StateMemberQz = 9,
	StateMemberAccBiasX = 10,
	StateMemberAccBiasY = 11,
	StateMemberAccBiasZ = 12,
	StateMemberGyroBiasX = 13,
	StateMemberGyroBiasY = 14,
	StateMemberGyroBiasZ = 15
};

constexpr int NOMINAL_STATE_SIZE{16};
constexpr int NOMINAL_POSITION_STATE_OFFSET = StateMemberX;
constexpr int NOMINAL_VELOCITY_STATE_OFFSET = StateMemberVx;
constexpr int NOMINAL_QUATERNION_STATE_OFFSET = StateMemberQw;
constexpr int NOMINAL_ACC_BIAS_STATE_OFFSET = StateMemberAccBiasX;
constexpr int NOMINAL_GYRO_BIAS_STATE_OFFSET = StateMemberGyroBiasX;
constexpr int NOMINAL_POSITION_STATE_SIZE{3};
constexpr int NOMINAL_VELOCITY_STATE_SIZE{3};
constexpr int NOMINAL_QUATERNION_STATE_SIZE{4};
constexpr int NOMINAL_ACC_BIAS_SIZE{3};
constexpr int NOMINAL_GYRO_BIAS_SIZE{3};

enum ErrorStateMembers  //(15x1)
{
	ErrorStateMemberX = 0,
	ErrorStateMemberY = 1,
	ErrorStateMemberZ = 2,
	ErrorStateMemberVx = 3,
	ErrorStateMemberVy = 4,
	ErrorStateMemberVz = 5,
	ErrorStateMemberDeltaRoll = 6,
	ErrorStateMemberDeltaPitch = 7,
	ErrorStateMemberDeltaYaw = 8,
	ErrorStateMemberAccBiasX = 9,
	ErrorStateMemberAccBiasY = 10,
	ErrorStateMemberAccBiasZ = 11,
	ErrorStateMemberGyroBiasX = 12,
	ErrorStateMemberGyroBiasY = 13,
	ErrorStateMemberGyroBiasZ = 14
};

constexpr double GRAVITY{ 9.80665 };
constexpr int ERROR_STATE_SIZE{15};
constexpr int ERROR_POSITION_STATE_OFFSET = ErrorStateMemberX;
constexpr int ERROR_VELOCITY_STATE_OFFSET = ErrorStateMemberVx;
constexpr int ERROR_DANGLE_STATE_OFFSET = ErrorStateMemberDeltaRoll;
constexpr int ERROR_ACC_BIAS_STATE_OFFSET = ErrorStateMemberAccBiasX;
constexpr int ERROR_GYRO_BIAS_STATE_OFFSET = ErrorStateMemberGyroBiasX;
constexpr int ERROR_POSITION_STATE_SIZE{3};
constexpr int ERROR_VELOCITY_STATE_SIZE{3};
constexpr int ERROR_DANGLE_STATE_SIZE{3};
constexpr int ERROR_ACC_BIAS_SIZE{3};
constexpr int ERROR_GYRO_BIAS_SIZE{3};




struct AdandGQGD {
	Matrix<double, ERROR_STATE_SIZE, ERROR_STATE_SIZE> Ad; //(15x15)
	Matrix<double, ERROR_STATE_SIZE, ERROR_STATE_SIZE> GQGD; //(15x15)
};

struct StatePredictions {
	Matrix<double,NOMINAL_STATE_SIZE,1> xNominalPrediction;				 // (16x1)
	Matrix<double,ERROR_STATE_SIZE,ERROR_STATE_SIZE> pPrediction;		 // (15x15)
};

struct InjectionStates {
	Matrix<double, NOMINAL_STATE_SIZE,1> xInject;				  // (16x1) 
	Matrix<double, ERROR_STATE_SIZE,ERROR_STATE_SIZE> pInject;  // (15x15)
};

struct InnovationPressureStates {
	//Matrix<double,1,1> pressureInnovation;				// (1x1)
	double pressureInnovation;
	Matrix<double,1,1> pressureInnovationCovariance;	// (1x1)
	Matrix<double,1,ERROR_STATE_SIZE> pressureH;		// (1x15)
};

struct InnovationDVLStates {
	MatrixXd DVLInnovation;
	MatrixXd DVLInnovationCovariance;
	MatrixXd DVLH;
};



class ESKF
{
public:

	explicit ESKF(const Matrix3d& Racc, const Matrix3d& RaccBias, const Matrix3d& Rgyro, const Matrix3d& RgyroBias,const double& pgyroBias,const double& paccBias, const Matrix3d& Sa, const Matrix3d& Sg,const Matrix3d& Sdvl, const Matrix3d& Sinc);

	VectorXd predictNominal(const VectorXd& xnominal,const Vector3d& accRectifiedMeasurements,const Vector3d& gyroRectifiedmeasurements,const double& Ts);
	MatrixXd Aerr(const VectorXd& xnominal, const Vector3d& accRectifiedMeasurements, const Vector3d& gyroRectifiedmeasurements);
	MatrixXd Gerr(const VectorXd& xnominal);
	AdandGQGD discreteErrorMatrix(const VectorXd& xnominal,const Vector3d& accRectifiedMeasurements,const Vector3d& gyroRectifiedmeasurements, const double& Ts);
	MatrixXd predictCovariance(const VectorXd& xnominal,const MatrixXd& P,const Vector3d& accRectifiedMeasurements,const Vector3d& gyroRectifiedmeasurements,const double& Ts);
	StatePredictions predict(const VectorXd& xnominal,const MatrixXd& P, Vector3d zAccMeasurements, Vector3d zGyroMeasurements,const double& Ts);
	InjectionStates inject(const VectorXd& xnominal,const VectorXd& deltaX,const MatrixXd& P);

	// DVL
	InnovationDVLStates innovationDVL(const VectorXd& xnominal,const MatrixXd& P,const Vector3d& zDVLvel,const Matrix3d& RDVL);
	InjectionStates updateDVL(const VectorXd& xnominal,const MatrixXd& P,const Vector3d& zDVLvel,const Matrix3d& RDVL);

	// Pressure sensor
	InnovationPressureStates innovationPressureZ(const VectorXd& xnominal,const MatrixXd& P,const double& zPressureZpos,const MatrixXd& RpressureZ);
	InjectionStates updatePressureZ(const VectorXd& xnominal,const MatrixXd& P,const double& zPressureZpos, const MatrixXd& RpressureZ);

private:
	double pgyroBias;
	double paccBias;

	const Vector3d gravity{ 0,0,GRAVITY };

	
	Matrix3d Racc;		// Acceleration measurements covariance (3x3)
	Matrix3d RaccBias;  // Acceleration bias driving noise covariance (3x3)
	Matrix3d Rgyro;		// Gyro measurements covariance (3x3)
	Matrix3d RgyroBias; // Gyro bias driving noise covariance (3x3)

	MatrixXd D; // Diagonal block matrix with measurement covariances


	// Correction matricies
	Matrix3d Sa; // Accelerometer
	Matrix3d Sg; // Gyro
	Matrix3d Sdvl; // DVL
	Matrix3d Sinc; // Inclinometer
	double SpressureZ; // Pressure
};

