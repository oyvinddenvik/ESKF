#include "ESKF.h"

#include <utility>
#include "common.h"

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
  , Sa_{ eulerToRotationMatrix(parameters.S_a) }
  , Sg_{ eulerToRotationMatrix(parameters.S_g) }
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
  std::cout << "S_a: " << parameters.S_a << std::endl;
  std::cout << "S_g: " << parameters.S_g << std::endl;
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

  MatrixXd A_err_discretized(ERROR_STATE_SIZE, ERROR_STATE_SIZE);
  MatrixXd A_err(ERROR_STATE_SIZE, ERROR_STATE_SIZE);
  MatrixXd identityMatrix(ERROR_STATE_SIZE, ERROR_STATE_SIZE);
  Matrix3d R_wdeltaT = Matrix3d::Zero();
  A_err.setZero();
  identityMatrix.setIdentity();

  R_wdeltaT = AngularErrorMatrix(gyroRectifiedmeasurements, Ts);

  A_err = Aerr(xnominal, accRectifiedMeasurements, gyroRectifiedmeasurements);

  A_err_discretized = identityMatrix + (A_err * Ts) + (0.5 * A_err * A_err * Ts * Ts) +
                      ((1.0 / 6.0) * A_err * A_err * A_err * Ts * Ts * Ts);

  A_err_discretized.block<3, 3>(6, 6) = R_wdeltaT.transpose();

  // Execution time
  // auto end = std::chrono::steady_clock::now();

  // auto diff = end - start;

  // std::cout <<"Aerr_discretized: " <<std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

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
  errorMatrix.Ad = AerrDiscretizedSecondOrder(xnominal, accRectifiedMeasurements, gyroRectifiedmeasurements, Ts);
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

InnovationParameters ESKF::innovationPressureZ(const VectorXd& xnominal, const MatrixXd& P, const double& zPressureZpos,
                                               const MatrixXd& RpressureZ)
{
  double zValue{ 1 };
  // Measurement Matrix
  MatrixXd Hx{ (MatrixXd(1, NOMINAL_STATE_SIZE) << MatrixXd::Zero(1, 2), zValue, MatrixXd::Zero(1, 16)).finished() };

  MatrixXd Q_deltaT{ (MatrixXd(4, 3) << -xnominal(7), -xnominal(8), -xnominal(9), xnominal(6), -xnominal(9),
                      xnominal(8), xnominal(9), xnominal(6), -xnominal(7), -xnominal(8), xnominal(7), xnominal(6))
                       .finished() };
  Q_deltaT *= 0.5;

  MatrixXd X_deltaX{ MatrixXd::Identity(NOMINAL_STATE_SIZE, ERROR_STATE_SIZE) };
  X_deltaX.block<4, 3>(6, 6) = Q_deltaT;

  InnovationParameters pressureStates(1);
  pressureStates.jacobianOfErrorStates = Hx * X_deltaX;
  pressureStates.measurementStates(0) = zPressureZpos - xnominal(StateMemberZ);
  pressureStates.measurementCovariance =
    (pressureStates.jacobianOfErrorStates * P * pressureStates.jacobianOfErrorStates.transpose()) + RpressureZ;

  return pressureStates;
}

void ESKF::updatePressureZ(const double& zPressureZpos, const MatrixXd& RpressureZ)
{
  InnovationParameters pressureStates =
    innovationPressureZ(optimizationParameters_.X, optimizationParameters_.P, zPressureZpos, RpressureZ);

  // ESKF Update step
  MatrixXd kalmanGain = optimizationParameters_.P * pressureStates.jacobianOfErrorStates.transpose() *
                        pressureStates.measurementCovariance.inverse();
  MatrixXd deltaX = kalmanGain * pressureStates.measurementStates;
  MatrixXd pUpdate =
    (identityMatrix_ - (kalmanGain * pressureStates.jacobianOfErrorStates)) * optimizationParameters_.P;
  optimizationParameters_ = inject(optimizationParameters_.X, deltaX, pUpdate);
}

InnovationParameters ESKF::innovationDVL(const VectorXd& xnominal, const MatrixXd& P, const Vector3d& zDVLvel,
                                         const Matrix3d& RDVL) const
{
  Vector3d vel_world = xnominal.block<NOMINAL_VELOCITY_STATE_SIZE, 1>(NOMINAL_VELOCITY_STATE_OFFSET, 0);
  Quaterniond nominalQuaternion{ optimizationParameters_.X(StateMemberQw), optimizationParameters_.X(StateMemberQx),
                                 optimizationParameters_.X(StateMemberQy), optimizationParameters_.X(StateMemberQz) };

  Matrix3d Hv{ nominalQuaternion.conjugate() };  // R_world_to_body

  // Todo: use analyic version from Sola paper
  MatrixXd jacobianMatrix = jacobianFdOfDVL(Hv * vel_world, nominalQuaternion.conjugate(), 0.000000001, vel_world);

  MatrixXd Hx{
    (MatrixXd(3, NOMINAL_STATE_SIZE) << Matrix3d::Zero(), Hv, jacobianMatrix, MatrixXd::Zero(3, 9)).finished()
  };

  MatrixXd Q_deltaT{ (MatrixXd(4, 3) << -nominalQuaternion.x(), -nominalQuaternion.y(), -nominalQuaternion.z(),
                      nominalQuaternion.w(), -nominalQuaternion.z(), nominalQuaternion.y(), nominalQuaternion.z(),
                      nominalQuaternion.w(), -nominalQuaternion.x(), -nominalQuaternion.y(), nominalQuaternion.x(),
                      nominalQuaternion.w())
                       .finished() };
  Q_deltaT *= 0.5;

  MatrixXd X_deltaX{ MatrixXd::Identity(NOMINAL_STATE_SIZE, ERROR_STATE_SIZE) };
  X_deltaX.block<4, 3>(6, 6) = Q_deltaT;

  InnovationParameters dvlStates(3);
  dvlStates.jacobianOfErrorStates = Hx * X_deltaX;
  dvlStates.measurementStates = zDVLvel - Hv * vel_world;
  dvlStates.measurementCovariance =
    (dvlStates.jacobianOfErrorStates * P * dvlStates.jacobianOfErrorStates.transpose()) + RDVL;

  return dvlStates;
}

void ESKF::updateDVL(const Vector3d& zDVLvel, const Matrix3d& RDVL)
{
  InnovationParameters DVLstates{ 3 };
  DVLstates = innovationDVL(optimizationParameters_.X, optimizationParameters_.P, zDVLvel, RDVL);

  MatrixXd kalmanGain =
    optimizationParameters_.P * DVLstates.jacobianOfErrorStates.transpose() * DVLstates.measurementCovariance.inverse();
  MatrixXd deltaX = kalmanGain * DVLstates.measurementStates;
  MatrixXd pUpdate = (identityMatrix_ - (kalmanGain * DVLstates.jacobianOfErrorStates)) * optimizationParameters_.P;

  optimizationParameters_ = inject(optimizationParameters_.X, deltaX, pUpdate);
}
