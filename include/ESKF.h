#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;

constexpr double GRAVITY{ 9.80665 };
constexpr int ROWSOFPOSITION{ 3 };
constexpr int ROWSOFVELOCITY{ 3 };
constexpr int ROWSOFQUATERNIONS{ 4 };
constexpr int ROWSOFACCBIAS{ 3 };
constexpr int ROWSOFGYROBIAS{ 3 };
constexpr int NumberOfNominalStates{ 16 };
constexpr int NumberOfErrorStates{ 15 };

struct AdandGQGD {
	Matrix<double, 15, 15> Ad;
	Matrix<double, 15, 15> GQGD;
};

struct StatePredictions {
	Matrix<double,NumberOfNominalStates,1> xNominalPrediction;				 // (16x1)
	Matrix<double,NumberOfErrorStates,NumberOfErrorStates> pPrediction;		 // (15x15)
};

struct InjectionStates {
	Matrix<double, NumberOfNominalStates,1> xInject;				  // (16x1) 
	Matrix<double, NumberOfErrorStates,NumberOfErrorStates> pInject;  // (15x15)
};

struct InnovationPressureStates {
	//Matrix<double,1,1> pressureInnovation;				// (1x1)
	double pressureInnovation;
	Matrix<double,1,1> pressureInnovationCovariance;	// (1x1)
	Matrix<double,1,NumberOfErrorStates> pressureH;		// (1x15)
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

