#include"common.h"
#include"ESKF.h"

using namespace Eigen;


int main()
{
    Matrix3d A;
    A << 1, 2, 3, 4, 5, 6, 7, 8, 9;
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

    std::cout << D << std::endl;


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