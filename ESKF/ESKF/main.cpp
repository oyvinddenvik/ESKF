#include"common.h"
#include"ESKF.h"
#include <Eigen/Core>

using namespace Eigen;


int main()
{
    int i = 3;

    std::cout << sqrt(i) << std::endl;
    Vector4d quatProduct;
    Vector4d quatRight;
    VectorXd quatLeft(4);
    quatLeft << 32, 1, 10, 11;
    quatRight << 0.5, 3, 6 , 10;
    VectorXd quatTest(4);
    VectorXd quatTest3(4);
    quatTest << 10, 20, 30, 40;
    quatTest3 << 5, 6, 7, 8;
    quatProduct.setZero();
    MatrixXd rotationMatrix;

    //std::cout << quatLeft.sum() << std::endl;
    
    std::cout << quatLeft.array().square().sum() << std::endl;
    //std::cout << quatLeft.block<3,1>(13,0) << std::endl;

    rotationMatrix = quaternion2Rotationmatrix(quatLeft);
    std::cout << rotationMatrix << std::endl;

    //quatProduct = quaternionHamiltonProduct(quatLeft, quatRight);

    //std::cout << quatProduct << std::endl;

    //quatProduct << quatLeft(0) * quatRight(0) - quatLeft.segment(1, 3).transpose() * quatRight.segment(1, 3),
      //  quatRight(0)* quatLeft.segment(1, 3) + quatLeft(0) * quatRight.segment(1, 3) + crossProductMatrix(quatLeft.segment(1, 3)) * quatRight.segment(1, 3);


    //Vector3d quatTest2;
    //quatTest.setZero();
    //quatTest << 1, 2, 3;
    //quatProduct.setZero();

    //quatProduct << 0, quatTest;

    //quatProduct.segment(1, 3) = quatRight;

    //quatProduct(1) = 40;

    //std::cout << quatProduct << std::endl << std::endl;
    //std::cout << quatProduct.segment(1, 3) << std::endl;
    //quatTest2 = quatProduct.segment(1, 3);
    //std::cout << quatProduct(1) << std::endl;
    //std::cout << quatTest2 << std::endl;

    //quatProduct.setZero();
    //quatProduct.count();

    //std::cout << quatProduct << std::endl;
    //std::cout << quatProduct.size() << std::endl;


    Matrix3d A;
    A << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    A.setIdentity();

    //std::cout << A << std::endl;

    DiagonalMatrix<double, 3> B{ 4,5,6 };
    DiagonalMatrix<double, 3> C{ 7,8,9 };
    DiagonalMatrix<double, 3> E{ 10,11,12 };

    MatrixXd D;

    D = blk3x3Diag(A, B, C, E);

    /*
    D.block<3, 3>(0, 0) = A;
    D.block<3, 3>(3, 3) = B;
    D.block<3, 3>(6, 6) = C;
    D.block<3, 3>(9, 9) = E;
    */

    //std::cout << D << std::endl;


    /*

    //const Vector3d gravity{1,1,GRAVITY};

    //std::cout << gravity << std::endl;

    Matrix3d b;
    Vector3d a{ 1, 2, 3 };
    std::cout << a << std::endl;
    b = crossProductMatrix(a);
    std::cout << b << std::endl;
    */



    /*
 
    MatrixXd a(3, 3);
    MatrixXd b(3, 3);
    MatrixXd c;

    c = blkdiag(a, 2);

    std::cout << c << std::endl;

    */
 
}