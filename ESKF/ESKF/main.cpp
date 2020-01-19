#include"common.h"
#include"ESKF.h"
#include <Eigen/Core>


using namespace Eigen;

struct structTest {
    Matrix<double, 3, 3> Ad;
    Matrix<double, 3, 3> Bd;
};

structTest returnMultipleMatricies()
{
    structTest testMatrix;

    testMatrix.Ad.setZero();
    testMatrix.Bd.setZero();

    testMatrix.Ad << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    testMatrix.Bd << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    return testMatrix;
}




int main()
{
    double Ts{ 0.0080 };
    Vector3d accRectifiedMeasurements = Vector3d::Zero();
    Vector3d gyroRectifiedMeasurements = Vector3d::Zero();
    VectorXd xnominal(16);
    Matrix3d Racc = Matrix3d::Zero();
    Matrix3d RaccBias = Matrix3d::Zero();
    Matrix3d Rgyro = Matrix3d::Zero();
    Matrix3d RgyroBias = Matrix3d::Zero();
    double pgyroBias{0.0001};
    double paccBias{ 0.0001 };
    Matrix3d Sa = Matrix3d::Zero();
    Matrix3d Sg = Matrix3d::Zero();
    Matrix3d Sdvl = Matrix3d::Zero();
    Matrix3d Sinc = Matrix3d::Zero();
    double SpressureZ;

    xnominal << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

    accRectifiedMeasurements << 1, 2, 3;
    gyroRectifiedMeasurements << 1, 2, 3;

    Racc << 3, 0, 0,
            0, 3, 0,
            0, 0, 3;

    RaccBias << 6e-5, 0, 0,
                0,6e-5 , 0,
                0, 0,6e-5;

    Rgyro << 0.0120, 0, 0,
            0, 0.0120, 0,
            0, 0, 0.0120;

    RgyroBias << 3e-7, 0, 0,
        0, 3e-7, 0,
        0, 0, 3e-7;


    ESKF eskf{ Racc,RaccBias,Rgyro,RgyroBias,pgyroBias,paccBias,Sa,Sg,Sdvl,Sinc};

    std::cout << eskf.predictNominal(xnominal,accRectifiedMeasurements,gyroRectifiedMeasurements,Ts)<<std::endl;
  

 
}