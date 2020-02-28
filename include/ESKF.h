#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <unsupported/Eigen/MatrixFunctions>
#include "common.h"



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
	StateMemberGyroBiasZ = 15,
	StateMemberGravityX = 16,
	StateMemberGravityY = 17,
	StateMemberGravityZ = 18
};

constexpr int NOMINAL_STATE_SIZE{19};  //16
constexpr int ERROR_STATE_SIZE{ 18 };
constexpr int NOMINAL_POSITION_STATE_OFFSET = StateMemberX;
constexpr int NOMINAL_VELOCITY_STATE_OFFSET = StateMemberVx;
constexpr int NOMINAL_QUATERNION_STATE_OFFSET = StateMemberQw;
constexpr int NOMINAL_ACC_BIAS_STATE_OFFSET = StateMemberAccBiasX;
constexpr int NOMINAL_GYRO_BIAS_STATE_OFFSET = StateMemberGyroBiasX;
constexpr int NOMINAL_GRAVITY_STATE_OFFSET = StateMemberGravityX;
constexpr int NOMINAL_POSITION_STATE_SIZE{3};
constexpr int NOMINAL_VELOCITY_STATE_SIZE{3};
constexpr int NOMINAL_QUATERNION_STATE_SIZE{4};
constexpr int NOMINAL_ACC_BIAS_SIZE{3};
constexpr int NOMINAL_GYRO_BIAS_SIZE{3};
constexpr int NOMINAL_GRAVITY_SIZE{3};


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
	Matrix<double,NOMINAL_VELOCITY_STATE_SIZE,1> DVLInnovation;
	Matrix<double,NOMINAL_VELOCITY_STATE_SIZE,NOMINAL_VELOCITY_STATE_SIZE > DVLInnovationCovariance;
	Matrix<double,NOMINAL_VELOCITY_STATE_SIZE,ERROR_STATE_SIZE> DVLH;
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
	Vector3d S_a;
	Vector3d S_g;
	//Matrix<double,3,3> S_a;
	//Matrix<double,3,3> S_g;
	Matrix<double,3,3> S_dvl;
	Matrix<double,3,3> S_inc;
	Matrix<double,NOMINAL_STATE_SIZE,1> initial_pose;
	Matrix<double,ERROR_STATE_SIZE,ERROR_STATE_SIZE> initial_covariance;
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
	AdandGQGD discreteErrorMatrix(const VectorXd& xnominal,const Vector3d& accRectifiedMeasurements,const Vector3d& gyroRectifiedmeasurements,const double& Ts, const Matrix3d& Racc, const Matrix3d& Rgyro);
	MatrixXd predictCovariance(const VectorXd& xnominal,const MatrixXd& P,const Vector3d& accRectifiedMeasurements,const Vector3d& gyroRectifiedmeasurements,const double& Ts, const Matrix3d& Racc, const Matrix3d& Rgyro);
	void predict(Vector3d zAccMeasurements, Vector3d zGyroMeasurements,const double& Ts, const Matrix3d& Racc, const Matrix3d& Rgyro);
	StatesAndErrorCovariance inject(const VectorXd& xnominal,const VectorXd& deltaX,const MatrixXd& P);

	// Lower computation time implementation
	MatrixXd AerrDiscretizedFirstOrder(const VectorXd& xnominal,const Vector3d& accRectifiedMeasurements,const Vector3d& gyroRectifiedmeasurements,const double& Ts);
	MatrixXd AerrDiscretizedSecondOrder(const VectorXd& xnominal,const Vector3d& accRectifiedMeasurements,const Vector3d& gyroRectifiedmeasurements,const double& Ts);
	MatrixXd AerrDiscretizedThirdOrder(const VectorXd& xnominal,const Vector3d& accRectifiedMeasurements,const Vector3d& gyroRectifiedmeasurements,const double& Ts);
	MatrixXd Fi();
	Matrix3d AngularErrorMatrix(const Vector3d& gyroRectifiedmeasurements,const double& Ts);

	// DVL
	InnovationDVLStates innovationDVL(const VectorXd& xnominal,const MatrixXd& P,const Vector3d& zDVLvel,const Matrix3d& RDVL);
	void updateDVL(const Vector3d& zDVLvel,const Matrix3d& RDVL);

	// Pressure sensor
	InnovationPressureStates innovationPressureZ(const VectorXd& xnominal,const MatrixXd& P,const double& zPressureZpos,const MatrixXd& RpressureZ);
	void updatePressureZ(const double& zPressureZpos, const MatrixXd& RpressureZ);


	void setParametersInESKF(const parametersInESKF& parameters);

	
	const inline Vector4d getQuaternion() const
	{
		Vector4d quaternion_NED = Vector4d::Zero();
		//Vector4d quaternion_NED_To_ENU = Vector4d::Zero();
		Vector4d quaternion_ENU = Vector4d::Zero();

		
		//quaternion_NED_To_ENU <<0, -0.70711, -0.70711, -0;

		quaternion_NED = poseStates_.block<NOMINAL_QUATERNION_STATE_SIZE, 1>(NOMINAL_QUATERNION_STATE_OFFSET, 0);

		//quaternion_ENU = quaternionHamiltonProduct(quaternion_NED_To_ENU, quaternion_NED);

		quaternion_ENU << quaternion_NED(0),       // w
						  quaternion_NED(2),       // y
						  quaternion_NED(1),       // x
						  -1.0*quaternion_NED(3);  // -z
		//std::cout<<quaternion_NED_To_ENU.conjugate()<<std::endl;

		if(use_ENU_)
		{
			return quaternion_ENU;
		}
		else
		{
			return quaternion_NED;
		}
		
	}

	const inline Vector3d getPosition() const
	{
		Matrix3d R_ned_to_enu = Matrix3d::Zero();
		Vector3d position = Vector3d::Zero();

		R_ned_to_enu << 0,1,0,
						1,0,0,
						0,0,-1;
		 
		position = poseStates_.block<NOMINAL_POSITION_STATE_SIZE,1>(NOMINAL_POSITION_STATE_OFFSET,0);

		if(use_ENU_)
		{
			return R_ned_to_enu*position;
		}
		else
		{
			return position;
		}
		
	}
	const inline Vector3d getVelocity() const
	{
		Matrix3d R_ned_to_enu = Matrix3d::Zero();
		Vector3d velocity = Vector3d::Zero();

		R_ned_to_enu << 0,1,0,
						1,0,0,
						0,0,-1;
		 
		velocity = poseStates_.block<NOMINAL_VELOCITY_STATE_SIZE, 1>(NOMINAL_VELOCITY_STATE_OFFSET, 0);

		//return R_ned_to_enu*velocity;
		
		if(use_ENU_)
		{
			return R_ned_to_enu*velocity;
		}
		else
		{
			return velocity;
		}
		
		 
	}

	const inline Vector3d getGravity() const
	{
		return poseStates_.block<NOMINAL_GRAVITY_SIZE,1>(NOMINAL_GRAVITY_STATE_OFFSET,0); 
	}

	const inline VectorXd getPose() const
	{
		return poseStates_;
	}

	const inline MatrixXd getErrorCovariance() const
	{
		return errorStateCovariance_;
	}

	
	


private:
	bool use_ENU_;
	double pgyroBias_;
	double paccBias_;

	//const Vector3d gravity_{ 0,0,GRAVITY };

	
	Matrix3d Racc_;		// Acceleration measurements covariance (3x3)
	Matrix3d RaccBias_;  // Acceleration bias driving noise covariance (3x3)
	Matrix3d Rgyro_;		// Gyro measurements covariance (3x3)
	Matrix3d RgyroBias_; // Gyro bias driving noise covariance (3x3)

	MatrixXd D_; // Diagonal block matrix with measurement covariances


	// Correction matricies
	Matrix3d Sa_; // Accelerometer
	Matrix3d Sg_; // Gyro
	Matrix3d Sdvl_; // DVL
	Matrix3d Sinc_; // Inclinometer
	double SpressureZ_; // Pressure

	VectorXd poseStates_;
	MatrixXd errorStateCovariance_;


	// Execution time
	std::vector<double> execution_time_vector_; 
	bool publish_execution_time_;


};

