#pragma once

#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

struct AdandGQGD {
	MatrixXd Ad;
	MatrixXd GQGD;
};

struct StatePredictions {
	MatrixXd xNominalPrediction; // ()
	MatrixXd pPrediction;		 // ()
};

struct InjectionStates {
	VectorXd xInject;  // (16x1) 
	MatrixXd pInject;  // (15x15)
};

struct InnovationPressureStates {
	MatrixXd pressureInnovation;		   // (1x1)
	MatrixXd pressureInnovationCovariance; // (1x1)
	MatrixXd pressureH;					   // (1x15)
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
	ESKF(Matrix3d Racc, Matrix3d RaccBias, Matrix3d Rgyro, Matrix3d RgyroBias, double pgyroBias, double paccBias, Matrix3d Sa, Matrix3d Sg, Matrix3d Sdvl, Matrix3d Sinc, double SpressureZ);

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

	const Vector3d gravity = Vector3d::Constant(0, 0, 9.80665);

	
	Matrix3d Racc;		// Acceleration measurements covariance (3x3)
	Matrix3d RaccBias;  // Acceleration bias driving noise covariance (3x3)
	Matrix3d Rgyro;		// Gyro measurements covariance (3x3)
	Matrix3d RgyroBias; // Gyro bias driving noise covariance (3x3)

	MatrixXd D;		    


	// Correction matricies
	Matrix3d Sa; // Accelerometer
	Matrix3d Sg; // Gyro
	Matrix3d Sdvl; // DVL
	Matrix3d Sinc; // Inclinometer
	double SpressureZ; // Pressure

};


