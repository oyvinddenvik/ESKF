#pragma once

#include <iostream>
#include <Eigen/Dense>


using namespace Eigen;


class ESKF
{
public:
	ESKF(Matrix3d Racc, Matrix3d RaccBias, Matrix3d Rgyro, Matrix3d RgyroBias, double pgyroBias, double paccBias, Matrix3d Sa, Matrix3d Sg, Matrix3d Sdvl, Matrix3d Sinc);

	VectorXd predictNominal(VectorXd xnominal, Vector3d accmes, Vector3d gyromes, double Ts);


private:
	double pgyroBias;
	double paccBias;

	Vector3d gravity;

	
	Matrix3d Racc;
	Matrix3d RaccBias;
	Matrix3d Rgyro;
	Matrix3d RgyroBias;

	// 
	MatrixXd D;


	// Correction matricies
	Matrix3d Sa; // accelerometer
	Matrix3d Sg; // gyro
	Matrix3d Sdvl; // DVL
	Matrix3d Sinc; // Inclinometer
	double SpressureZ; // Pressure


};


