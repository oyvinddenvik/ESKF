#include "ESKF.h"
#include <utility>
#include "common.h"
#include <vector>

using namespace eskf;

ESKF::ESKF(Matrix3d Racc, Matrix3d RaccBias, Matrix3d Rgyro, Matrix3d RgyroBias, double pgyroBias, double paccBias,
           Matrix3d Sa, Matrix3d Sg, Matrix3d Sdvl, Matrix3d Sinc)
  : pgyroBias_{ pgyroBias }
  , paccBias_{ paccBias }
  , Racc_{ std::move(Racc) }
  , RaccBias_{ std::move(RaccBias) }
  , Rgyro_{ std::move(Rgyro) }
  , RgyroBias_{ std::move(RgyroBias) }
  , D_{ blk3x3Diag(Racc_, Rgyro_, RaccBias_, RgyroBias_) }
  , Sa_{ std::move(Sa) }
  , Sg_{ std::move(Sg) }
  , Sdvl_{ std::move(Sdvl) }
  , Sinc_{ std::move(Sinc) }
{
  Vector3d initialPosition{ -1.96, 0.061, 0 };
  Vector3d initialVelocity{ 0, 0, 0 };
  Vector4d initialQuat{ 0.9935, 0, 0, -0.1134 };
  Vector3d initialAccBias{ 0, 0, 0 };
  Vector3d initialGyroBias{ 0, 0, 0 };
  Vector3d initialGravity{ 0, 0, 9.80665 };

  optimizationParameters_.X << initialPosition, initialVelocity, initialQuat, initialAccBias, initialGyroBias,
    initialGravity;

  Vector3d initialPPos{ 1e-9 * 1e-9, 1e-9 * 1e-9, 1e-9 * 1e-9 };
  Vector3d initialPVel{ 1e-4 * 1e-4, 1e-4 * 1e-4, 1e-4 * 1e-4 };
  Vector3d initialPDq{ 12e-3 * 12e-3, 12e-3 * 12e-3, 12e-3 * 12e-3 };
  Vector3d initialPAccBias{ 12e-9 * 12e-9, 12e-9 * 12e-9, 12e-9 * 12e-9 };
  Vector3d initialPGyroBias{ 3e-18 * 3e-18, 3e-18 * 3e-18, 3e-18 * 3e-18 };
  Vector3d initialPGravity{ 3e-18 * 3e-18, 3e-18 * 3e-18, 3e-18 * 3e-18 };

  optimizationParameters_.P.block<3, 3>(0, 0) = initialPPos.asDiagonal();
  optimizationParameters_.P.block<3, 3>(3, 3) = initialPVel.asDiagonal();
  optimizationParameters_.P.block<3, 3>(6, 6) = initialPDq.asDiagonal();
  optimizationParameters_.P.block<3, 3>(9, 9) = initialPAccBias.asDiagonal();
  optimizationParameters_.P.block<3, 3>(12, 12) = initialPGyroBias.asDiagonal();
  optimizationParameters_.P.block<3, 3>(15, 15) = initialPGravity.asDiagonal();
}



ESKF::ESKF(const parametersInESKF& parameters)
  : Racc_{ parameters.R_acc }
  , Rgyro_{ parameters.R_gyro }
  , RaccBias_{ parameters.R_accBias }
  , RgyroBias_{ parameters.R_gyroBias }
  , Sdvl_{ parameters.S_dvl }
  , Sinc_{ parameters.S_inc }
  , Sa_{ eulerToRotationMatrix(parameters.Sr_to_ned_accelerometer + parameters.Sr_accelerometer_aligment)} //+ eulerToRotationMatrix(parameters.Sr_accelerometer_aligment)}
  , Sg_{ eulerToRotationMatrix(parameters.Sr_to_ned_gyro + parameters.Sr_gyro_aligment)} //+ eulerToRotationMatrix(parameters.Sr_gyro_aligment)}
  , pgyroBias_{ parameters.pgyroBias }
  , paccBias_{ parameters.paccBias }
  , use_ENU_{ parameters.use_ENU }
  , D_{ blk3x3Diag(Racc_, Rgyro_, RaccBias_, RgyroBias_) }
  , optimizationParameters_{ parameters.initial_pose, parameters.initial_covariance }
{
  std::cout << "R_acc: " << parameters.R_acc << std::endl;
  std::cout << "R_gyro: " << parameters.R_gyro << std::endl;
  std::cout << "R_gyroBias: " << parameters.R_gyroBias << std::endl;
  std::cout << "R_accBias: " << parameters.R_accBias << std::endl;
  std::cout << "S_DVL: " << parameters.S_dvl << std::endl;
  std::cout << "Sr_to_ned_accelerometer: " << parameters.Sr_to_ned_accelerometer << std::endl;
  std::cout << "Sr_to_ned_gyro: " << parameters.Sr_to_ned_gyro << std::endl;
  std::cout << "Sr_accelerometer_alignment: "<<parameters.Sr_accelerometer_aligment << std::endl;
  std::cout << "Sr_gyro_alignment: "<<parameters.Sr_gyro_aligment << std::endl; 
  std::cout << "S_a: "<<Sa_<<std::endl;
  std::cout << "Sr_to_ned_accelerometer + Sr_accelerometer_alignment: " << parameters.Sr_to_ned_accelerometer + parameters.Sr_accelerometer_aligment<< std::endl;

}

VectorXd ESKF::predictNominal(const VectorXd& xnominal, const Vector3d& accRectifiedMeasurements,
                              const Vector3d& gyroRectifiedmeasurements, const double& Ts) const
{
  // Extract states
  Vector3d position = xnominal.block<NOMINAL_POSITION_STATE_SIZE, 1>(NOMINAL_POSITION_STATE_OFFSET, 0);
  Vector3d velocity = xnominal.block<NOMINAL_VELOCITY_STATE_SIZE, 1>(NOMINAL_VELOCITY_STATE_OFFSET, 0);
  Quaterniond quaternion{ optimizationParameters_.X(StateMemberQw), optimizationParameters_.X(StateMemberQx),
                          optimizationParameters_.X(StateMemberQy), optimizationParameters_.X(StateMemberQz) };
  Vector3d accBias = xnominal.block<NOMINAL_ACC_BIAS_SIZE, 1>(NOMINAL_ACC_BIAS_STATE_OFFSET, 0);
  Vector3d gyroBias = xnominal.block<NOMINAL_GYRO_BIAS_SIZE, 1>(NOMINAL_GYRO_BIAS_STATE_OFFSET, 0);
  Vector3d gravity = xnominal.block<NOMINAL_GRAVITY_SIZE, 1>(NOMINAL_GRAVITY_STATE_OFFSET, 0);

  // Predictions
  // quaternion product is overloaded: Vector3d = Quaternion * Vector3d <-- Vector4d = Quaternion * Vector4d *
  // Quaternion*
  position += (Ts * velocity) + (0.5 * Ts * Ts * ((quaternion * accRectifiedMeasurements) + gravity));
  velocity += Ts * ((quaternion * accRectifiedMeasurements) + gravity);

  // q <-- q * q{(omega_m - omega_b)*dt} sola (260c)
  // where q{omega} = Exp(omega)         sola (101)

  // TODO: Use Sophus libriary: exp(Tangent dTheta)
  Vector3d dTheta = Ts * gyroRectifiedmeasurements;
  double dAngle = sqrt(dTheta.array().square().sum());
  Quaterniond dq{ cos(dAngle / 2.0), sin(dAngle / 2.0) * (dTheta(0) / dAngle), sin(dAngle / 2.0) * (dTheta(1) / dAngle),
                  sin(dAngle / 2.0) * (dTheta(2) / dAngle) };
  quaternion *= dq;
  // TODO: Check if correct (this is the nominal update)
  accBias = exp(-1.0 * paccBias_ * Ts) * accBias;
  gyroBias = exp(-1.0 * pgyroBias_ * Ts) * gyroBias;
  gravity = gravity;

  // Normalize quaternion
  quaternion.normalize();

  // Concatenate
  VectorXd xNextnominal{ VectorXd::Zero(NOMINAL_STATE_SIZE) };
  xNextnominal << position, velocity, quaternion.w(), quaternion.x(), quaternion.y(), quaternion.z(), accBias, gyroBias,
    gravity;

  return xNextnominal;
}

MatrixXd ESKF::Aerr(const VectorXd& xnominal, const Vector3d& accRectifiedMeasurements,
                    const Vector3d& gyroRectifiedmeasurements) const
{
  Quaterniond quaternion{ optimizationParameters_.X(StateMemberQw), optimizationParameters_.X(StateMemberQx),
                          optimizationParameters_.X(StateMemberQy), optimizationParameters_.X(StateMemberQz) };
  Matrix3d rotationMatrix{ quaternion };

  MatrixXd A{ MatrixXd::Zero(ERROR_STATE_SIZE, ERROR_STATE_SIZE) };
  A.block<3, 3>(0, 3) = Matrix3d::Identity();
  A.block<3, 3>(3, 6) = -1.0 * rotationMatrix * crossProductMatrix(accRectifiedMeasurements);  // TODO: use Sophus
  A.block<3, 3>(3, 9) = -1.0 * rotationMatrix;
  A.block<3, 3>(6, 6) = -1.0 * crossProductMatrix(gyroRectifiedmeasurements);
  A.block<3, 3>(6, 12) = -Matrix3d::Identity();
  A.block<3, 3>(9, 9) = -1.0 * paccBias_ * Matrix3d::Identity();
  A.block<3, 3>(12, 12) = -1.0 * pgyroBias_ * Matrix3d::Identity();

  // Gravity estimation
  A.block<3, 3>(3, 15) = Matrix3d::Identity();

  // Bias corrections
  A.block<3, 3>(3, 9) = A.block<3, 3>(3, 9) * Sa_;
  A.block<3, 3>(6, 12) = A.block<3, 3>(6, 12) * Sg_;

  return A;
}

Matrix3d ESKF::AngularErrorMatrix(const Vector3d& gyroRectifiedmeasurements, const double& Ts) const
{
  return (
    Matrix3d::Identity() - (crossProductMatrix(gyroRectifiedmeasurements) * sin(Ts)) +
    (crossProductMatrix(gyroRectifiedmeasurements) * crossProductMatrix(gyroRectifiedmeasurements) * (1 - cos(Ts))));
}

MatrixXd ESKF::AerrDiscretizedFirstOrder(const VectorXd& xnominal, const Vector3d& accRectifiedMeasurements,
                                         const Vector3d& gyroRectifiedmeasurements, const double& Ts) const
{
  // auto start = std::chrono::steady_clock::now();
  // Compute first order approximation of matrix exponentional taylor expansion
  // I_15x15 +A_err*Ts

  MatrixXd A_err_discretized(ERROR_STATE_SIZE, ERROR_STATE_SIZE);
  MatrixXd A_err(ERROR_STATE_SIZE, ERROR_STATE_SIZE);
  MatrixXd identityMatrix(ERROR_STATE_SIZE, ERROR_STATE_SIZE);
  Matrix3d R_wdeltaT = Matrix3d::Zero();
  A_err.setZero();
  identityMatrix.setIdentity();

  R_wdeltaT = AngularErrorMatrix(gyroRectifiedmeasurements, Ts);

  A_err = Aerr(xnominal, accRectifiedMeasurements, gyroRectifiedmeasurements);

  A_err_discretized = identityMatrix + A_err * Ts;

  A_err_discretized.block<3, 3>(6, 6) = R_wdeltaT.transpose();

  // Execution time
  // auto end = std::chrono::steady_clock::now();

  // auto diff = end - start;

  // std::cout <<"Aerr_discretized: " <<std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

  return A_err_discretized;
}

MatrixXd ESKF::AerrDiscretizedSecondOrder(const VectorXd& xnominal, const Vector3d& accRectifiedMeasurements,
                                          const Vector3d& gyroRectifiedmeasurements, const double& Ts) const
{
  // Compute second order approximation of matrix exponentional taylor expansion
  // I_15x15 +A_err*Ts

  Matrix3d R_wdeltaT = AngularErrorMatrix(gyroRectifiedmeasurements, Ts);

  MatrixXd A_err = Aerr(xnominal, accRectifiedMeasurements, gyroRectifiedmeasurements);

  MatrixXd A_err_discretized = identityMatrix_ + A_err * Ts + 0.5 * A_err * A_err * Ts * Ts;

  A_err_discretized.block<3, 3>(6, 6) = R_wdeltaT.transpose();

  return A_err_discretized;
}

MatrixXd ESKF::AerrDiscretizedThirdOrder(const VectorXd& xnominal, const Vector3d& accRectifiedMeasurements,
                                         const Vector3d& gyroRectifiedmeasurements, const double& Ts) const
{
  // auto start = std::chrono::steady_clock::now();
  // Compute second order approximation of matrix exponentional taylor expansion
  // I_15x15 +A_err*Ts

  Matrix3d R_wdeltaT = AngularErrorMatrix(gyroRectifiedmeasurements, Ts);

  MatrixXd A_err = Aerr(xnominal, accRectifiedMeasurements, gyroRectifiedmeasurements);

  MatrixXd A_err_discretized = identityMatrix_ + (A_err * Ts) + (0.5 * A_err * A_err * Ts * Ts) +
                      ((1.0 / 6.0) * A_err * A_err * A_err * Ts * Ts * Ts);

  A_err_discretized.block<3, 3>(6, 6) = R_wdeltaT.transpose();

  return A_err_discretized;
}

MatrixXd ESKF::Fi()
{
  MatrixXd F_i(ERROR_STATE_SIZE, 12);
  Matrix3d identity_matrix_3x3 = Matrix3d::Zero();
  identity_matrix_3x3.setIdentity();
  F_i.setZero();

  F_i.block<3, 3>(3, 0) = identity_matrix_3x3;
  F_i.block<3, 3>(6, 3) = identity_matrix_3x3;
  F_i.block<3, 3>(9, 6) = identity_matrix_3x3;
  F_i.block<3, 3>(12, 9) = identity_matrix_3x3;

  return F_i;
}

MatrixXd ESKF::Gerr(const VectorXd& xnominal)
{
  Quaterniond quaternion{ xnominal(StateMemberQw), xnominal(StateMemberQx), xnominal(StateMemberQy),
                          xnominal(StateMemberQz) };

  MatrixXd Gerror{ MatrixXd::Zero(ERROR_STATE_SIZE, 12) };
  Gerror.block<3, 3>(3, 0) = -1.0 * quaternion.toRotationMatrix();
  Gerror.block<3, 3>(6, 3) = -1.0 * Matrix3d::Identity();
  Gerror.block<3, 3>(9, 6) = Matrix3d::Identity();
  Gerror.block<3, 3>(12, 9) = Matrix3d::Identity();

  return Gerror;
}

AdandGQGD ESKF::discreteErrorMatrix(const VectorXd& xnominal, const Vector3d& accRectifiedMeasurements,
                                    const Vector3d& gyroRectifiedmeasurements, const double& Ts, const Matrix3d& Racc,
                                    const Matrix3d& Rgyro) const
{
  // auto start = std::chrono::steady_clock::now();

  // Initilize
  AdandGQGD errorMatrix;
  // MatrixXd vanLoan(ERROR_STATE_SIZE*2, ERROR_STATE_SIZE*2);
  // MatrixXd vanLoanExponentional(ERROR_STATE_SIZE*2, ERROR_STATE_SIZE*2);
  // MatrixXd zeros(ERROR_STATE_SIZE, ERROR_STATE_SIZE);
  // MatrixXd A(ERROR_STATE_SIZE, ERROR_STATE_SIZE);
  // MatrixXd G(ERROR_STATE_SIZE, 12);
  // MatrixXd FDF(ERROR_STATE_SIZE,ERROR_STATE_SIZE);

  // A.setZero();
  // G.setZero();
  // vanLoanExponentional.setZero();
  // zeros.setZero();
  // vanLoan.setZero();
  // FDF.setZero();

  // D_ = blk3x3Diag(Racc, Rgyro, RaccBias_, RgyroBias_);

  // A = Aerr(xnominal, accRectifiedMeasurements, gyroRectifiedmeasurements);
  // G = Gerr(xnominal);

  // Calculate Van Loan
  // vanLoan << -1.0 * A, G* D_* G.transpose(),
  //			zeros, A.transpose();

  // vanLoan = vanLoan * Ts;

  // Computation time very slow
  // vanLoanExponentional= vanLoan.exp();

  // errorMatrix.Ad = vanLoanExponentional.block<ERROR_STATE_SIZE, ERROR_STATE_SIZE>(ERROR_STATE_SIZE,
  // ERROR_STATE_SIZE).transpose(); errorMatrix.GQGD = vanLoanExponentional.block<ERROR_STATE_SIZE,
  // ERROR_STATE_SIZE>(ERROR_STATE_SIZE, ERROR_STATE_SIZE).transpose() * vanLoanExponentional.block<ERROR_STATE_SIZE,
  // ERROR_STATE_SIZE>(0, ERROR_STATE_SIZE);
  errorMatrix.Ad = AerrDiscretizedThirdOrder(xnominal, accRectifiedMeasurements, gyroRectifiedmeasurements, Ts);
  errorMatrix.GQGD = F_i_ * D_ * F_i_.transpose();

  // Execution time
  // auto end = std::chrono::steady_clock::now();

  // auto diff = end - start;

  // std::cout <<"DiscreteErrorMatrix: " <<std::chrono::duration <double, std::milli> (diff).count() << " ms" <<
  // std::endl;

  return std::move(errorMatrix);
}

MatrixXd ESKF::predictCovariance(const VectorXd& xnominal, const MatrixXd& P, const Vector3d& accRectifiedMeasurements,
                                 const Vector3d& gyroRectifiedmeasurements, const double& Ts, const Matrix3d& Racc,
                                 const Matrix3d& Rgyro) const
{
  AdandGQGD errorMatrix =
    discreteErrorMatrix(xnominal, accRectifiedMeasurements, gyroRectifiedmeasurements, Ts, Racc, Rgyro);
  MatrixXd Pprediction = (errorMatrix.Ad * P * errorMatrix.Ad.transpose()) + errorMatrix.GQGD;

  return Pprediction;
}

void ESKF::predict(const Vector3d& zAccMeasurements, const Vector3d& zGyroMeasurements, const double& Ts,
                   const Matrix3d& Racc, const Matrix3d& Rgyro)
{
  

  Vector3d accBias = Sa_ * optimizationParameters_.X.block<NOMINAL_ACC_BIAS_SIZE, 1>(NOMINAL_ACC_BIAS_STATE_OFFSET, 0);
  Vector3d gyroBias =
    Sg_ * optimizationParameters_.X.block<NOMINAL_GYRO_BIAS_SIZE, 1>(NOMINAL_GYRO_BIAS_STATE_OFFSET, 0);

  Vector3d accelerationRectified = Sa_ * zAccMeasurements - accBias;
  Vector3d gyroRectified = Sg_ * zGyroMeasurements - gyroBias;

  optimizationParameters_.X = predictNominal(optimizationParameters_.X, accelerationRectified, gyroRectified, Ts);
  optimizationParameters_.P = predictCovariance(optimizationParameters_.X, optimizationParameters_.P,
                                                accelerationRectified, gyroRectified, Ts, Racc, Rgyro);

  
                                              
}

StatesAndErrorCovariance ESKF::inject(const VectorXd& xnominal, const VectorXd& deltaX, const MatrixXd& P) const
{
  Vector3d positionInjections =
    xnominal.block<NOMINAL_POSITION_STATE_SIZE, 1>(NOMINAL_POSITION_STATE_OFFSET, 0) + deltaX.block<3, 1>(0, 0);
  Vector3d velocityInjections =
    xnominal.block<NOMINAL_VELOCITY_STATE_SIZE, 1>(NOMINAL_VELOCITY_STATE_OFFSET, 0) + deltaX.block<3, 1>(3, 0);
  Quaterniond quatRight{ 1, deltaX(6) / 2.0, deltaX(7) / 2.0, deltaX(8) / 2.0 };
  Quaterniond nominalQuaternion{ xnominal(StateMemberQw), xnominal(StateMemberQx), xnominal(StateMemberQy),
                                 xnominal(StateMemberQz) };
  Quaterniond quaternionInjections = nominalQuaternion * quatRight;
  Vector3d accelerationBiasInjections =
    xnominal.block<NOMINAL_ACC_BIAS_SIZE, 1>(NOMINAL_ACC_BIAS_STATE_OFFSET, 0) + deltaX.block<3, 1>(9, 0);
  Vector3d gyroBiasInjections =
    xnominal.block<NOMINAL_GYRO_BIAS_SIZE, 1>(NOMINAL_GYRO_BIAS_STATE_OFFSET, 0) + deltaX.block<3, 1>(12, 0);

  // Gravity estimation
  Vector3d gravityInjections =
    xnominal.block<NOMINAL_GRAVITY_SIZE, 1>(NOMINAL_GRAVITY_STATE_OFFSET, 0) + deltaX.block<3, 1>(15, 0);

  // Normalize quaternion
  quaternionInjections.normalize();

  StatesAndErrorCovariance injections;
  injections.X << positionInjections, velocityInjections, quaternionInjections.w(), quaternionInjections.x(),
    quaternionInjections.y(), quaternionInjections.z(), accelerationBiasInjections, gyroBiasInjections,
    gravityInjections;

  MatrixXd Ginject{ MatrixXd::Identity(ERROR_STATE_SIZE, ERROR_STATE_SIZE) };
  Ginject.block<3, 3>(6, 6) = Matrix3d::Identity() - crossProductMatrix(0.5 * deltaX.block<3, 1>(6, 0));
  injections.P = Ginject * P * Ginject.transpose();

  return injections;
}

InnovationParameters ESKF::innovationPosition(const VectorXd& xnominal, const MatrixXd& P, const VectorXd& positionMeasurement, const MatrixXd& RPosition, const std::vector<int>& indexXYZ)
{
  assert(positionMeasurement.rows() <= 3 && "position Measurements rows are greater than 3");
  assert(indexXYZ.size() <= 3 && "indexXYZ size is greater than 3");

  // Measurement Matrix
  MatrixXd Hx{ (MatrixXd(positionMeasurement.rows(), NOMINAL_STATE_SIZE) << MatrixXd::Zero(positionMeasurement.rows(), NOMINAL_POSITION_STATE_SIZE-positionMeasurement.rows()), MatrixXd::Identity(positionMeasurement.rows(),positionMeasurement.rows()), MatrixXd::Zero(positionMeasurement.rows(), NOMINAL_STATE_SIZE-positionMeasurement.rows())).finished() };

  MatrixXd Q_deltaT{ (MatrixXd(4, 3) << -xnominal(7), -xnominal(8), -xnominal(9), xnominal(6), -xnominal(9),
                      xnominal(8), xnominal(9), xnominal(6), -xnominal(7), -xnominal(8), xnominal(7), xnominal(6))
                       .finished() };
  Q_deltaT *= 0.5;

  MatrixXd X_deltaX{ MatrixXd::Identity(NOMINAL_STATE_SIZE, ERROR_STATE_SIZE) };
  X_deltaX.block<4, 3>(6, 6) = Q_deltaT;

  InnovationParameters positionStates(positionMeasurement.rows());
  positionStates.jacobianOfErrorStates = Hx * X_deltaX;
  for(u_int8_t i{0}; i<positionMeasurement.rows();++i)
  {
    positionStates.measurementStates(i) = positionMeasurement(i) - xnominal(i);
  }
  positionStates.measurementCovariance =
    (positionStates.jacobianOfErrorStates * P * positionStates.jacobianOfErrorStates.transpose()) + RPosition;

  return positionStates;
}

void ESKF::updatePosition(const VectorXd& rawPosition, const MatrixXd& RPosition)
{
  assert(rawPosition.rows() <= 3 && "position Measurements rows are greater than 3");

  InnovationParameters positionStates =
    innovationPosition(optimizationParameters_.X, optimizationParameters_.P, rawPosition, RPosition);

  // ESKF Update step
  MatrixXd kalmanGain = optimizationParameters_.P * positionStates.jacobianOfErrorStates.transpose() *
                        positionStates.measurementCovariance.inverse();
  MatrixXd deltaX = kalmanGain * positionStates.measurementStates;
  MatrixXd pUpdate =
    (identityMatrix_ - (kalmanGain * positionStates.jacobianOfErrorStates)) * optimizationParameters_.P;
  optimizationParameters_ = inject(optimizationParameters_.X, deltaX, pUpdate);
}

InnovationParameters ESKF::innovationVelocity(const VectorXd& xnominal, const MatrixXd& P, const VectorXd& velocityMeasurement,
                                         const MatrixXd& RVel, const std::vector<int>& indexXYZ) const
{
  
  assert(velocityMeasurement.rows() <= 3 && "velocity Measurements rows are greater than 3");
  assert(indexXYZ.size() <= 3 && "indexXYZ size is greater than 3");

  Matrix<double,3,1> vel_world;
  vel_world.setZero();

  Quaterniond nominalQuaternion{ optimizationParameters_.X(StateMemberQw), optimizationParameters_.X(StateMemberQx),
                                 optimizationParameters_.X(StateMemberQy), optimizationParameters_.X(StateMemberQz) };

  Matrix3d Hv{ nominalQuaternion.conjugate() };  // R_world_to_body

  // TODO: use analyic version from Sola paper
  MatrixXd jacobianMatrix = jacobianFdOfDVL(Hv * vel_world, nominalQuaternion.conjugate(), 0.000000001, vel_world);


  std::vector<MatrixXd> getColumnOfHv;
  std::vector<MatrixXd> getRowsOfJacobianMatrix;
  for (auto i : indexXYZ)
  {
    vel_world(i) = xnominal(NOMINAL_VELOCITY_STATE_OFFSET + i);

    getColumnOfHv.push_back(Hv.block(0,i,3,1));
    getRowsOfJacobianMatrix.push_back(jacobianMatrix.block(i,0,1,4));

    //MatrixXd getColumnOfHv{Hv.block(0,0,3,i)};
    //MatrixXd getRowsOfJacobianMatrix{jacobianMatrix.block(0,0,i,4)};
  };

  MatrixXd columnsOfHv;

  for(u_int8_t i{0}; i<getColumnOfHv.size();++i)
  {
    columnsOfHv.block<3,1>(0,i) = getColumnOfHv[i];
  }

  MatrixXd rowsOfJacobianMatrix;

  for(u_int8_t i{0}; i<getRowsOfJacobianMatrix.size();++i)
  {
    rowsOfJacobianMatrix.block<1,3>(i,0) = getRowsOfJacobianMatrix[i];
  }

  MatrixXd Hx{
  (MatrixXd(indexXYZ.size(), NOMINAL_STATE_SIZE) << MatrixXd::Zero(velocityMeasurement.rows(), NOMINAL_POSITION_STATE_SIZE-indexXYZ.size()), columnsOfHv.transpose(), rowsOfJacobianMatrix, MatrixXd::Zero(indexXYZ.size(), NOMINAL_STATE_SIZE-indexXYZ.size())).finished()
  };

  /*
  for(u_int8_t i{0}; i<velocityMeasurement.rows();++i)
  {
    vel_world(i) = xnominal(NOMINAL_VELOCITY_STATE_OFFSET + i);
  }
  
  //xnominal.block<NOMINAL_VELOCITY_STATE_SIZE, 1>(NOMINAL_VELOCITY_STATE_OFFSET, 0);
  MatrixXd Hx{
    (MatrixXd(velocityMeasurement.rows(), NOMINAL_STATE_SIZE) << MatrixXd::Zero(velocityMeasurement.rows(), NOMINAL_POSITION_STATE_SIZE-velocityMeasurement.rows()), getColumnOfHv.transpose(), getRowsOfJacobianMatrix, MatrixXd::Zero(velocityMeasurement.rows(), NOMINAL_STATE_SIZE-velocityMeasurement.rows())).finished()
  };
  */

  MatrixXd Q_deltaT{ (MatrixXd(4, 3) << -nominalQuaternion.x(), -nominalQuaternion.y(), -nominalQuaternion.z(),
                      nominalQuaternion.w(), -nominalQuaternion.z(), nominalQuaternion.y(), nominalQuaternion.z(),
                      nominalQuaternion.w(), -nominalQuaternion.x(), -nominalQuaternion.y(), nominalQuaternion.x(),
                      nominalQuaternion.w())
                       .finished() };
  Q_deltaT *= 0.5;

  MatrixXd X_deltaX{ MatrixXd::Identity(NOMINAL_STATE_SIZE, ERROR_STATE_SIZE) };
  X_deltaX.block<4, 3>(6, 6) = Q_deltaT;

  InnovationParameters velocityStates(velocityMeasurement.rows());
  velocityStates.jacobianOfErrorStates = Hx * X_deltaX;

  //velocityStates.measurementStates = velocityMeasurement - Hv * vel_world;


  velocityStates.measurementStates.block(0,0,velocityMeasurement.rows(),1) = velocityMeasurement - Hv.block(0,0,velocityMeasurement.rows(),3) * vel_world.block(0,0,velocityMeasurement.rows(),1);
  

  for (auto i : indexXYZ)
  {
    velocityStates.measurementStates.block(0,0,indexXYZ.size(),1) = velocityMeasurement - columnsOfHv * vel_world.block(i,0,indexXYZ.size(),1);
  }

  velocityStates.measurementCovariance =
    (velocityStates.jacobianOfErrorStates * P * velocityStates.jacobianOfErrorStates.transpose()) + RVel;

  return velocityStates;
}

void ESKF::updateVelocity(const VectorXd& velocityMeasurement, const MatrixXd& RVel, const std::vector<int>& indexXYZ)
{
  InnovationParameters velocitystates{ velocityMeasurement.rows() };
  velocitystates = innovationVelocity(optimizationParameters_.X, optimizationParameters_.P, velocityMeasurement, RVel, indexXYZ);

  MatrixXd kalmanGain =
    optimizationParameters_.P * velocitystates.jacobianOfErrorStates.transpose() * velocitystates.measurementCovariance.inverse();
  MatrixXd deltaX = kalmanGain * velocitystates.measurementStates;
  MatrixXd pUpdate = (identityMatrix_ - (kalmanGain * velocitystates.jacobianOfErrorStates)) * optimizationParameters_.P;

  optimizationParameters_ = inject(optimizationParameters_.X, deltaX, pUpdate);
}


void updatePose(const VectorXd& rawMeasurement, const MatrixXd& RPose,const Matrix<bool,3,3> sensorStateConfig)
{

  // Get position from sensor_state_config


  std::vector<int> a;

  for(u_int8_t i{0}; i<sensorStateConfig.rows();++i)
  {
      for(u_int8_t j{0}; i < sensorStateConfig.cols(); ++j){
        if sensorStateConfig(i,j)
        {
          a.push_back(i*3+j);
        }   
      }
  }








}
