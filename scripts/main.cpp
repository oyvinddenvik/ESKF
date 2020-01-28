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
    InjectionStates injectionDVLTesting;
    InnovationDVLStates dvlstatetesting;
    Vector3d ZdvlValues = Vector3d::Zero();
    Matrix3d RDVL = Matrix3d::Zero();
    InjectionStates injectionPressureTesting;
    StatePredictions statetesting;
    InnovationPressureStates pressureTesting;
    VectorXd deltaX(15);
    InjectionStates injectiontesting;
    MatrixXd P(15, 15);
    AdandGQGD testing;
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
    double pressureValue{ 0.5 };
    MatrixXd RpressureZ(1, 1);

    P.setIdentity();

    xnominal << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
    deltaX << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;

    accRectifiedMeasurements << 1, 2, 3;
    gyroRectifiedMeasurements << 1, 2, 3;

    ZdvlValues << 1, 2, 3;
    RDVL.setIdentity();

    RDVL = RDVL * 2.0;


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

    RpressureZ << 0.5;

    Sa.setIdentity();
    Sg.setIdentity();

    ESKF eskf{ Racc,RaccBias,Rgyro,RgyroBias,pgyroBias,paccBias,Sa,Sg,Sdvl,Sinc};
    
    Racc = Racc * Ts;

    pressureTesting = eskf.innovationPressureZ(xnominal, P, pressureValue, RpressureZ);
    injectionPressureTesting = eskf.updatePressureZ(xnominal, P, pressureValue, RpressureZ);
    dvlstatetesting = eskf.innovationDVL(xnominal, P, ZdvlValues, RDVL);
    injectionDVLTesting = eskf.updateDVL(xnominal, P, ZdvlValues, RDVL);

    //std::cout << injectionDVLTesting.pInject << std::endl;
    std::cout << injectionDVLTesting.xInject << std::endl;

    //std::cout << dvlstatetesting.DVLH << std::endl;
    //std::cout << dvlstatetesting.DVLInnovation << std::endl;
    //std::cout << dvlstatetesting.DVLInnovationCovariance << std::endl;



    //std::cout << pressureTesting.pressureH << std::endl;
    //std::cout << pressureTesting.pressureInnovation << std::endl;
    //std::cout << pressureTesting.pressureInnovationCovariance << std::endl;

    //std::cout << injectionPressureTesting.xInject << std::endl;
    //std::cout << injectionPressureTesting.pInject << std::endl;


 

    //statetesting = eskf.predict(xnominal, P, accRectifiedMeasurements, gyroRectifiedMeasurements, Ts);
    //injectiontesting = eskf.inject(xnominal, deltaX, P);
    //testing = eskf.discreteErrorMatrix(xnominal, accRectifiedMeasurements, gyroRectifiedMeasurements, Ts);

    //std::cout << eskf.predictNominal(xnominal,accRectifiedMeasurements,gyroRectifiedMeasurements,Ts)<<std::endl;
    //std::cout << eskf.Aerr(xnominal,accRectifiedMeasurements,gyroRectifiedMeasurements)<<std::endl;
    //std::cout << eskf.Gerr(xnominal)<<std::endl;
    //std::cout << testing.Ad <<std::endl;
    //std::cout << testing.GQGD << std::endl;
    //eskf.discreteErrorMatrix(xnominal, accRectifiedMeasurements, gyroRectifiedMeasurements, Ts);
    
    //std::cout << eskf.predictCovariance(xnominal,P,accRectifiedMeasurements,gyroRectifiedMeasurements,Ts) <<std::endl;
    //std::cout << injectiontesting.xInject << std::endl;
    //std::cout << statetesting.xNominalPrediction << std::endl;
 
}