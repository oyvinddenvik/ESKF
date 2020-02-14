#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <unsupported/Eigen/MatrixFunctions>

//#include "ros_node.h"

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
constexpr int ERROR_STATE_SIZE{ 15 };
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


constexpr double GRAVITY{ 9.80665 };
constexpr int DEFAULT_IMU_RATE{125};




struct AdandGQGD {
	Matrix<double, ERROR_STATE_SIZE,ERROR_STATE_SIZE> Ad;  //(15x15)
	Matrix<double, ERROR_STATE_SIZE,ERROR_STATE_SIZE> GQGD; //(15x15)
};


struct StatesAndErrorCovariance {
	Matrix<double,NOMINAL_STATE_SIZE,1> X;
	Matrix<double,ERROR_STATE_SIZE,ERROR_STATE_SIZE> P;
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

struct parametersInESKF 
{
    Matrix<double,3,3> R_acc;
    Matrix<double,3,3> R_accBias;
	Matrix<double,3,3> R_gyro;
	Matrix<double,3,3> R_gyroBias;
	Matrix<double,3,3> R_dvl;
	Matrix<double,1,1> R_pressureZ;
	double pgyroBias;
	double paccBias;
	//Vector3d S_a;
	Matrix<double,3,3> S_a;
	Matrix<double,3,3> S_g;
	Matrix<double,3,3> S_dvl;
	Matrix<double,3,3> S_inc;
	
	bool use_ENU;
};



class ESKF
{
public:

	ESKF();

	explicit ESKF(const Matrix3d& Racc, const Matrix3d& RaccBias, const Matrix3d& Rgyro, const Matrix3d& RgyroBias,const double& pgyroBias,const double& paccBias, const Matrix3d& Sa, const Matrix3d& Sg,const Matrix3d& Sdvl, const Matrix3d& Sinc);

	VectorXd predictNominal(const VectorXd& xnominal,const Vector3d& accRectifiedMeasurements,const Vector3d& gyroRectifiedmeasurements,const double& Ts);
	MatrixXd Aerr(const VectorXd& xnominal, const Vector3d& accRectifiedMeasurements, const Vector3d& gyroRectifiedmeasurements);
	MatrixXd Gerr(const VectorXd& xnominal);
	AdandGQGD discreteErrorMatrix(const VectorXd& xnominal,const Vector3d& accRectifiedMeasurements,const Vector3d& gyroRectifiedmeasurements, const double& Ts);
	MatrixXd predictCovariance(const VectorXd& xnominal,const MatrixXd& P,const Vector3d& accRectifiedMeasurements,const Vector3d& gyroRectifiedmeasurements,const double& Ts);
	void predict(Vector3d zAccMeasurements, Vector3d zGyroMeasurements,const double& Ts);
	//StatePredictions predict(const VectorXd& xnominal,const MatrixXd& P, Vector3d zAccMeasurements, Vector3d zGyroMeasurements,const double& Ts);
	StatesAndErrorCovariance inject(const VectorXd& xnominal,const VectorXd& deltaX,const MatrixXd& P);

	// DVL
	InnovationDVLStates innovationDVL(const VectorXd& xnominal,const MatrixXd& P,const Vector3d& zDVLvel,const Matrix3d& RDVL);
	void updateDVL(const Vector3d& zDVLvel,const Matrix3d& RDVL);

	// Pressure sensor
	InnovationPressureStates innovationPressureZ(const VectorXd& xnominal,const MatrixXd& P,const double& zPressureZpos,const MatrixXd& RpressureZ);
	void updatePressureZ(const double& zPressureZpos, const MatrixXd& RpressureZ);


	void setParametersInESKF(const parametersInESKF& parameters);

	const inline Matrix3d getS_a() const
	{
		return Sa;
	}

	const inline Vector3d getPositionInENU() const
	{
		Matrix3d R_ned_to_enu = Matrix3d::Zero();
		Vector3d position = Vector3d::Zero();

		R_ned_to_enu << 0,1,0,
						1,0,0,
						0,0,-1;
		 
		position = poseStates.block<NOMINAL_POSITION_STATE_SIZE,1>(NOMINAL_POSITION_STATE_OFFSET,0);

		return R_ned_to_enu*position;
		/*
		if(use_ENU_)
		{
			return R_ned_to_enu*position;
		}
		else
		{
			return position;
		}
		*/
	}
	const inline Vector3d getVelocityInENU() const
	{
		Matrix3d R_ned_to_enu = Matrix3d::Zero();
		Vector3d velocity = Vector3d::Zero();

		R_ned_to_enu << 0,1,0,
						1,0,0,
						0,0,-1;
		 
		velocity = poseStates.block<NOMINAL_VELOCITY_STATE_SIZE, 1>(NOMINAL_VELOCITY_STATE_OFFSET, 0);

		return R_ned_to_enu*velocity;
		/*
		if(use_ENU_)
		{
			return R_ned_to_enu*velocity;
		}
		else
		{
			return velocity;
		}
		*/
		 
	}

	const inline VectorXd getPose() const
	{
		return poseStates;
	}

	const inline MatrixXd getErrorCovariance() const
	{
		return errorStateCovariance;
	}

	
	


private:
	bool use_ENU_;
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

	VectorXd poseStates;
	MatrixXd errorStateCovariance;

};

