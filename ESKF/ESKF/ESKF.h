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

	// IMU and ESKF implementation

	//ESKF();
	ESKF(Matrix3d Racc, Matrix3d RaccBias, Matrix3d Rgyro, Matrix3d RgyroBias, double pgyroBias, double paccBias, Matrix3d Sa, Matrix3d Sg, Matrix3d Sdvl, Matrix3d Sinc);

	VectorXd predictNominal(VectorXd xnominal, Vector3d accRectifiedMeasurements, Vector3d gyroRectifiedmeasurements, double Ts);
	MatrixXd Aerr(VectorXd xnominal, Vector3d accRectifiedMeasurements, Vector3d gyroRectifiedmeasurements);
	MatrixXd Gerr(VectorXd xnominal);
	AdandGQGD discreteErrorMatrix(VectorXd xnominal, Vector3d accRectifiedMeasurements, Vector3d gyroRectifiedmeasurements, double Ts);
	MatrixXd predictCovariance(VectorXd xnominal, MatrixXd P, Vector3d accRectifiedMeasurements, Vector3d gyroRectifiedmeasurements, double Ts);
	StatePredictions predict(VectorXd xnominal, MatrixXd P, VectorXd zAccMeasurements, VectorXd zGyroMeasurements, double Ts);
	InjectionStates inject(VectorXd xnominal, VectorXd deltaX, MatrixXd P);


	// DVL
	InnovationDVLStates innovationDVL(VectorXd xnominal, MatrixXd P, Vector3d zDVLvel, MatrixXd RDVL);
	InjectionStates updateDVL(VectorXd xnominal, MatrixXd P, Vector3d zDVLvel, MatrixXd RDVL);


	// Pressure sensor
	InnovationPressureStates innovationPressureZ(VectorXd xnominal, MatrixXd P, double zPressureZpos, MatrixXd RpressureZ);
	InjectionStates updatePressureZ(VectorXd xnominal, MatrixXd P, double zPressureZpos, MatrixXd RpressureZ);

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

