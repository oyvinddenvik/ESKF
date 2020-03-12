//#include"common.h"
#include "ESKF.h"
#include <Eigen/Core>
#include "ros_node.h"
#include "chrono"

using namespace eskf;
using namespace Eigen;


ESKF_Node::ESKF_Node(const ros::NodeHandle& nh, const ros::NodeHandle& pnh)
  : nh_{ pnh }
  , init_{ false }
  , eskf_{ loadParametersFromYamlFile() }
  , publish_execution_time_{true}
  //,eskf_{R_ACC,R_ACCBIAS,R_GYRO,R_GYROBIAS,P_GYRO_BIAS,P_ACC_BIAS,returnStaticRotationFromIMUtoBodyFrame(roll_pitch_yaw_NED_and_alignment_corrected),returnStaticRotationFromIMUtoBodyFrame(roll_pitch_yaw_NED_and_alignment_corrected),S_DVL,S_INC}
{
  R_dvl_.setZero();
  R_pressureZ_.setZero();
  std::string imu_topic{ "" };
  std::string dvl_topic{ "" };
  std::string pressureZ_topic{ "" };
  int publish_rate{ 125 };

  // const parametersInESKF parameters = loadParametersFromYamlFile();
  // eskf_.setParametersInESKF(parameters);

  setRdvlFromYamlFile(R_dvl_);
  setRpressureZFromYamlFile(R_pressureZ_);
  setIMUTopicNameFromYaml(imu_topic);
  setDVLTopicNameFromYawl(dvl_topic);
  setPressureZTopicNameFromYaml(pressureZ_topic);
  setPublishrateFromYaml(publish_rate);

  ROS_INFO("Subscribing to IMU topic: %s", imu_topic.c_str());
  // Subscribe to IMU
  subscribeIMU_ = nh_.subscribe<sensor_msgs::Imu>(imu_topic, 1000, &ESKF_Node::imuCallback, this,
                                                  ros::TransportHints().tcpNoDelay(true));
  // Subscribe to DVL
  ROS_INFO("Subscribing to DVL: %s", dvl_topic.c_str());
  subcribeDVL_ = nh_.subscribe<nav_msgs::Odometry>(dvl_topic, 1000, &ESKF_Node::dvlCallback, this,
                                                   ros::TransportHints().tcpNoDelay(true));
  // Subscribe to Pressure sensor
  ROS_INFO("Subscribing to pressure sensor: %s", pressureZ_topic.c_str());
  subscribePressureZ_ = nh_.subscribe<nav_msgs::Odometry>(pressureZ_topic, 1000, &ESKF_Node::pressureZCallback, this,
                                                          ros::TransportHints().tcpNoDelay(true));

  ROS_INFO("Publishing State");
  publishPose_ = nh_.advertise<nav_msgs::Odometry>("pose", 1);

  pubTImer_ = nh_.createTimer(ros::Duration(1.0f / publish_rate), &ESKF_Node::publishPoseState, this);
}

// IMU Subscriber
void ESKF_Node::imuCallback(const sensor_msgs::Imu::ConstPtr& imu_Message_data)
{

  auto start = std::chrono::steady_clock::now();

  double Ts{ 0 };
  int imu_publish_rate{ DEFAULT_IMU_RATE };
  Vector3d raw_acceleration_measurements = Vector3d::Zero();
  Vector3d raw_gyro_measurements = Vector3d::Zero();
  Matrix3d R_acc = Matrix3d::Zero();
  Matrix3d R_gyro = Matrix3d::Zero();

  Ts = (1.0 / imu_publish_rate);
  raw_acceleration_measurements << imu_Message_data->linear_acceleration.x, imu_Message_data->linear_acceleration.y,
    imu_Message_data->linear_acceleration.z;

  for (size_t i = 0; i < R_acc.rows(); i++)
  {
    for (size_t j = 0; j < R_acc.cols(); j++)
    {
      R_acc(i, j) = imu_Message_data->linear_acceleration_covariance[3 * i + j];
    }
  }

  raw_gyro_measurements << imu_Message_data->angular_velocity.x, imu_Message_data->angular_velocity.y,
    imu_Message_data->angular_velocity.z;

  for (size_t i = 0; i < R_gyro.rows(); i++)
  {
    for (size_t j = 0; j < R_gyro.cols(); j++)
    {
      R_gyro(i, j) = imu_Message_data->angular_velocity_covariance[3 * i + j];
    }
  }

  // ROS_INFO("Acceleration_x: %f",imu_Message_data->linear_acceleration.x);

  

  eskf_.predict(raw_acceleration_measurements, raw_gyro_measurements, Ts, R_acc, R_gyro);


  // Execution time
	auto end = std::chrono::steady_clock::now();

	auto diff = end - start;

	auto diff_in_ms = std::chrono::duration <double, std::milli> (diff).count();

  //std::cout<<diff_in_ms<<std::endl;

	
	if(execution_time_vector_.size() == 1000 && publish_execution_time_ == true)
	{
		std::cout<<"Max value: "<<maxOfVector(execution_time_vector_)<<std::endl;
		std::cout<<"Mean: "<<meanOfVector(execution_time_vector_)<<std::endl;
		std::cout<<"STD: "<<stanardDeviationOfVector(execution_time_vector_)<<std::endl;
		publish_execution_time_ = false;
	}
	else
	{
		execution_time_vector_.push_back(diff_in_ms);
	}
}

// DVL subscriber
void ESKF_Node::dvlCallback(const nav_msgs::Odometry::ConstPtr& dvl_Message_data)
{
  Vector3d raw_dvl_measurements = Vector3d::Zero();
  Matrix3d R_dvl = Matrix3d::Zero();

  raw_dvl_measurements << dvl_Message_data->twist.twist.linear.x, dvl_Message_data->twist.twist.linear.y,
    dvl_Message_data->twist.twist.linear.z;

  for (size_t i = 0; i < R_dvl.rows(); i++)
  {
    for (size_t j = 0; j < R_dvl.cols(); j++)
    {
      R_dvl(i, j) = dvl_Message_data->twist.covariance[3 * i + j];
    }
  }

  // ROS_INFO("Velocity_z: %f",dvl_Message_data->twist.twist.linear.z);
  eskf_.updateDVL(raw_dvl_measurements, R_dvl_);
}

// PressureZ subscriber
void ESKF_Node::pressureZCallback(const nav_msgs::Odometry::ConstPtr& pressureZ_Message_data)
{
  Matrix<double, 1, 1> RpressureZ;
  const double raw_pressure_z = pressureZ_Message_data->pose.pose.position.z;

  RpressureZ(0) = pressureZ_Message_data->pose.covariance[0];

  // std::cout<<RpressureZ<<std::endl;
  // const double R_pressureZ = 2.2500;

  eskf_.updatePressureZ(raw_pressure_z, R_pressureZ_);
}

void ESKF_Node::publishPoseState(const ros::TimerEvent&)
{
  nav_msgs::Odometry odom_msg;
  static size_t trace_id{ 0 };

  const Vector3d& position = eskf_.getPosition();
  const Vector3d& velocity = eskf_.getVelocity();
  Quaterniond quaternion = eskf_.getQuaternion();
  const Vector3d& gravity = eskf_.getGravity();

  // const VectorXd& pose = eskf_.getPose();
  // const MatrixXd& errorCovariance = eskf_.getErrorCovariance();

  odom_msg.header.frame_id = "/eskf_link";
  odom_msg.header.seq = trace_id++;
  odom_msg.header.stamp = ros::Time::now();
  odom_msg.pose.pose.position.x = position(StateMemberX);  // NEDpose(StateMemberX);
  odom_msg.pose.pose.position.y = position(StateMemberY);  // NEDpose(StateMemberY)
  odom_msg.pose.pose.position.z = position(StateMemberZ);  // NEDpose(StateMemberZ);
  odom_msg.twist.twist.linear.x = velocity(0);             // NEDpose(StateMemberVx);
  odom_msg.twist.twist.linear.y = velocity(1);             // NEDpose(StateMemberVy);
  odom_msg.twist.twist.linear.z = velocity(2);             // NEDpose(StateMemberVz);
  odom_msg.pose.pose.orientation.w = quaternion.w();       // pose(StateMemberQw);
  odom_msg.pose.pose.orientation.x = quaternion.x();       // pose(StateMemberQx);
  odom_msg.pose.pose.orientation.y = quaternion.y();       // pose(StateMemberQy);
  odom_msg.pose.pose.orientation.z = quaternion.z();       // pose(StateMemberQz);
  // odom_msg.pose.covariance

  // loadParametersFromYamlFile();

  /*
  Quat q;
  EulerAngles euler;

  q.w = odom_msg.pose.pose.orientation.w;
  q.x = odom_msg.pose.pose.orientation.x;
  q.y = odom_msg.pose.pose.orientation.y;
  q.z = odom_msg.pose.pose.orientation.z;

  euler = fromQuaternionToEulerAngles(q);

  */
  //auto euler = quaternion.toRotationMatrix().eulerAngles(0, 1, 2) * 180.0 / 3.14;
  //std::cout << "Roll pitch yaw" << euler << std::endl;
  //std::cout << "gravity_X: " << gravity(0) << std::endl;
  //std::cout << "gravity_Y: " << gravity(1) << std::endl;
  //std::cout << "gravity_Z: " << gravity(2) << std::endl;

  /*
  // Position covariance
  for(size_t i = 0; i < NOMINAL_POSITION_STATE_SIZE;i++)
  {
      for(size_t j = 0; j<NOMINAL_POSITION_STATE_SIZE;j++)
      {
          odom_msg.pose.covariance[NOMINAL_POSITION_STATE_SIZE*i+j] = errorCovariance(i,j);
      }
  }

  for(size_t i = 0; i < NOMINAL_VELOCITY_STATE_SIZE; i++)
  {
      for(size_t j = 0; j< NOMINAL_VELOCITY_STATE_SIZE; j++)
      {
          odom_msg.pose.covariance[(NOMINAL_VELOCITY_STATE_SIZE*i+j)+8] = errorCovariance(i +
  NOMINAL_VELOCITY_STATE_OFFSET, j + NOMINAL_VELOCITY_STATE_OFFSET);
      }
  }
  */

  /*
   for (size_t i = 0; i < POSE_SIZE; i++)
    {
      for (size_t j = 0; j < POSE_SIZE; j++)
      {
        message.pose.covariance[POSE_SIZE * i + j] = estimateErrorCovariance(i, j);
      }
    }
  */

  // ROS_INFO("StateX: %f",odom_msg.pose.pose.position.x);
  publishPose_.publish(odom_msg);

  // std::cout<<pose<<std::endl;
  // std::cout<<std::endl;
}



parametersInESKF ESKF_Node::loadParametersFromYamlFile()
{
  parametersInESKF parameters;
  parameters.R_acc.setZero();
  parameters.R_accBias.setZero();
  parameters.R_gyro.setZero();
  parameters.R_gyroBias.setZero();
  parameters.R_dvl.setZero();
  parameters.R_pressureZ.setZero();
  parameters.Sr_to_ned_accelerometer.setZero();
  parameters.Sr_to_ned_gyro.setZero();
  parameters.Sr_accelerometer_aligment.setZero();
  parameters.Sr_gyro_aligment.setZero();
  parameters.Sr_to_ned_dvl.setZero();
  parameters.Sr_dvl_alignment.setZero();
  parameters.S_inc.setZero();
  parameters.S_dvl.setZero();
  parameters.paccBias = 0;
  parameters.pgyroBias = 0;
  parameters.use_ENU = false;
  parameters.initial_pose.setZero();
  parameters.initial_covariance.setZero();
  // parameters.initial_pose(16);
  // parameters.initial_pose.setZero();

  // std::cout<<parameters.initial_pose<<std::endl;
  // Eigen::MatrixXd R_acc(3,3);
  // R_acc.setZero();

  XmlRpc::XmlRpcValue R_dvlConfig;
  XmlRpc::XmlRpcValue R_pressureZConfig;
  XmlRpc::XmlRpcValue R_accConfig;
  XmlRpc::XmlRpcValue R_accBiasConfig;
  XmlRpc::XmlRpcValue R_gyroConfig;
  XmlRpc::XmlRpcValue R_gyroBiasConfig;
  XmlRpc::XmlRpcValue St_to_ned_accelerometer_Config;
  XmlRpc::XmlRpcValue St_to_ned_gyro_Config;
  XmlRpc::XmlRpcValue St_accelerometer_alignment_Config;
  XmlRpc::XmlRpcValue St_gyro_alignment_Config;
  XmlRpc::XmlRpcValue St_to_ned_DVL_Config;
  XmlRpc::XmlRpcValue St_dvl_alignment_Config;
  XmlRpc::XmlRpcValue S_gConfig;
  XmlRpc::XmlRpcValue S_dvlConfig;
  XmlRpc::XmlRpcValue S_incConfig;
  XmlRpc::XmlRpcValue initialPoseConfig;
  XmlRpc::XmlRpcValue initialCovarianceConfig;

  if (ros::param::has("/R_acc"))
  {
    ros::param::get("/R_acc", R_accConfig);
    int matrix_size = parameters.R_acc.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      for (int j = 0; j < matrix_size; j++)
      {
        std::ostringstream ostr;
        ostr << R_accConfig[matrix_size * i + j];
        std::istringstream istr(ostr.str());
        istr >> parameters.R_acc(i, j);
      }
    }
  }
  else
  {
    ROS_FATAL("No measurement covariance for accelerometer (R_acc) set in parameter file!");
    ROS_BREAK();
  }

  if (ros::param::has("/R_accBias"))
  {
    ros::param::get("/R_accBias", R_accBiasConfig);

    int matrix_size = parameters.R_accBias.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      for (int j = 0; j < matrix_size; j++)
      {
        std::ostringstream ostr;
        ostr << R_accBiasConfig[matrix_size * i + j];
        std::istringstream istr(ostr.str());
        istr >> parameters.R_accBias(i, j);
      }
    }
  }
  else
  {
    ROS_FATAL("No measurement covariance for accelerometer bias (R_accBias) set in parameter file!");
    ROS_BREAK();
  }

  if (ros::param::has("/R_gyro"))
  {
    ros::param::get("/R_gyro", R_gyroConfig);
    int matrix_size = parameters.R_gyro.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      for (int j = 0; j < matrix_size; j++)
      {
        std::ostringstream ostr;
        ostr << R_gyroConfig[matrix_size * i + j];
        std::istringstream istr(ostr.str());
        istr >> parameters.R_gyro(i, j);
      }
    }
  }
  else
  {
    ROS_FATAL("No measurement covariance for gyro (R_gyro) set in parameter file");
    ROS_BREAK();
  }

  if (ros::param::has("/R_gyroBias"))
  {
    ros::param::get("/R_gyroBias", R_gyroBiasConfig);
    int matrix_size = parameters.R_gyroBias.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      for (int j = 0; j < matrix_size; j++)
      {
        std::ostringstream ostr;
        ostr << R_gyroBiasConfig[matrix_size * i + j];
        std::istringstream istr(ostr.str());
        istr >> parameters.R_gyroBias(i, j);
      }
    }
  }
  else
  {
    ROS_FATAL("No measurement covariance for gyro bias (R_gyroBias) set in parameter file");
    ROS_BREAK();
  }

  if (ros::param::has("/sr_accelerometer_to_NED"))
  {
    ros::param::get("/sr_accelerometer_to_NED", St_to_ned_accelerometer_Config);
    int matrix_size = parameters.Sr_to_ned_accelerometer.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      std::ostringstream ostr;
      ostr <<St_to_ned_accelerometer_Config[i];
      std::istringstream istr(ostr.str());
      istr >> parameters.Sr_to_ned_accelerometer(i);
    }
  }
  else
  {
    ROS_FATAL("No static transform for accelerometer (St_acc) set in parameter file");
    ROS_BREAK();
  }

  if (ros::param::has("/sr_gyro_to_NED"))
  {
    ros::param::get("/sr_gyro_to_NED", St_to_ned_gyro_Config);
    int matrix_size = parameters.Sr_to_ned_gyro.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      std::ostringstream ostr;
      ostr << St_to_ned_gyro_Config[i];
      std::istringstream istr(ostr.str());
      istr >> parameters.Sr_to_ned_gyro(i);
    }
  }
  else
  {
    ROS_FATAL("No static transform for gyro (St_gyro) set in parameter file");
    ROS_BREAK();
  }

  if (ros::param::has("/sr_accelerometer_alignment"))
  {
    ros::param::get("/sr_accelerometer_alignment", St_accelerometer_alignment_Config);
    int matrix_size = parameters.Sr_accelerometer_aligment.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      std::ostringstream ostr;
      ostr << St_accelerometer_alignment_Config[i];
      std::istringstream istr(ostr.str());
      istr >> parameters.Sr_accelerometer_aligment(i);
    }
  }
  else
  {
    ROS_FATAL("No static transform for gyro (St_gyro) set in parameter file");
    ROS_BREAK();
  }



  if (ros::param::has("/sr_gyro_alignment"))
  {
    ros::param::get("/sr_gyro_alignment", St_gyro_alignment_Config);
    int matrix_size = parameters.Sr_gyro_aligment.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      std::ostringstream ostr;
      ostr << St_gyro_alignment_Config[i];
      std::istringstream istr(ostr.str());
      istr >> parameters.Sr_gyro_aligment(i);
    }
  }
  else
  {
    ROS_FATAL("No static transform for gyro (St_gyro) set in parameter file");
    ROS_BREAK();
  }

   if (ros::param::has("/sr_dvl_alignment"))
  {
    ros::param::get("/sr_dvl_alignment", St_dvl_alignment_Config);
    int matrix_size = parameters.Sr_dvl_alignment.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      std::ostringstream ostr;
      ostr << St_dvl_alignment_Config[i];
      std::istringstream istr(ostr.str());
      istr >> parameters.Sr_dvl_alignment(i);
    }
  }
  else
  {
    ROS_FATAL("No static rotation for alignment of Dvl set in parameter file");
    ROS_BREAK();
  }

   if (ros::param::has("/sr_dvl_to_NED"))
  {
    ros::param::get("/sr_dvl_to_NED", St_to_ned_DVL_Config);
    int matrix_size = parameters.Sr_to_ned_dvl.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      std::ostringstream ostr;
      ostr << St_to_ned_DVL_Config[i];
      std::istringstream istr(ostr.str());
      istr >> parameters.Sr_to_ned_dvl(i);
    }
  }
  else
  {
    ROS_FATAL("No static rotation for DVL (sr_dvl_to_NED) set in parameter file");
    ROS_BREAK();
  }




  if (ros::param::has("/St_dvl"))
  {
    ros::param::get("/St_dvl", S_dvlConfig);
    int matrix_size = parameters.S_dvl.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      for (int j = 0; j < matrix_size; j++)
      {
        std::ostringstream ostr;
        ostr << S_dvlConfig[matrix_size * i + j];
        std::istringstream istr(ostr.str());
        istr >> parameters.S_dvl(i, j);
      }
    }
  }
  else
  {
    ROS_FATAL("No static transform for DVL (St_dvl) set in parameter file");
    ROS_BREAK();
  }

  if (ros::param::has("/St_inc"))
  {
    ros::param::get("/St_inc", S_incConfig);
    int matrix_size = parameters.S_inc.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      for (int j = 0; j < matrix_size; j++)
      {
        std::ostringstream ostr;
        ostr << S_incConfig[matrix_size * i + j];
        std::istringstream istr(ostr.str());
        istr >> parameters.S_inc(i, j);
      }
    }
  }
  else
  {
    ROS_FATAL("No static transform for inclinometer (St_inc) set in parameter file");
    ROS_BREAK();
  }

  if (ros::param::has("/R_dvl"))
  {
    ros::param::get("/R_dvl", R_dvlConfig);
    int matrix_size = parameters.R_dvl.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      for (int j = 0; j < matrix_size; j++)
      {
        std::ostringstream ostr;
        ostr << R_dvlConfig[matrix_size * i + j];
        std::istringstream istr(ostr.str());
        istr >> parameters.R_dvl(i, j);
      }
    }
  }
  else
  {
    ROS_FATAL("No measurement covariance for DVL (R_dvl) set in parameter file");
    ROS_BREAK();
  }

  if (ros::param::has("/p_gyroBias"))
  {
    ros::param::get("/p_gyroBias", parameters.pgyroBias);
  }
  else
  {
    ROS_FATAL("No gyro bias driving noise (p_gyroBias) set in parameter file!");
    ROS_BREAK();
  }

  if (ros::param::has("/p_accBias"))
  {
    ros::param::get("/p_accBias", parameters.paccBias);
  }
  else
  {
    ROS_FATAL("No acceleration bias driving noise (p_accBias) set in parameter file!");
    ROS_BREAK();
  }

  if (ros::param::has("/publish_in_ENU"))
  {
    ros::param::get("/publish_in_ENU", parameters.use_ENU);
  }
  else
  {
    ROS_FATAL("No bool value set for publish_in_ENU in parameter file! ");
    ROS_BREAK();
  }

  if (ros::param::has("/R_pressureZ"))
  {
    ros::param::get("/R_pressureZ", R_pressureZConfig);
    int matrix_size = parameters.R_pressureZ.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      std::ostringstream ostr;
      ostr << R_pressureZConfig[i];
      std::istringstream istr(ostr.str());
      istr >> parameters.R_pressureZ(i);
    }
  }
  else
  {
    ROS_FATAL("No measurement covariance for pressure sensor (R_pressureZ) set in parameter file!");
    ROS_BREAK();
  }

  if (ros::param::has("/initial_pose"))
  {
    ros::param::get("/initial_pose", initialPoseConfig);
    int matrix_size = parameters.initial_pose.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      std::ostringstream ostr;
      ostr << initialPoseConfig[i];
      std::istringstream istr(ostr.str());
      istr >> parameters.initial_pose(i);
    }
  }
  else
  {
    ROS_FATAL("No initial pose (initial_pose) set in parameter file!");
    ROS_BREAK();
  }

  if (ros::param::has("/initial_covariance"))
  {
    ros::param::get("/initial_covariance", initialCovarianceConfig);
    int matrix_size = parameters.initial_covariance.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      for (int j = 0; j < matrix_size; j++)
      {
        std::ostringstream ostr;
        ostr << initialCovarianceConfig[matrix_size * i + j];
        std::istringstream istr(ostr.str());
        istr >> parameters.initial_covariance(i, j);
      }
    }
  }
  else
  {
    ROS_FATAL("No initial covariance (initial_covariance) set in parameter file!");
    ROS_BREAK();
  }

  return parameters;
}

void setIMUTopicNameFromYaml(std::string& imu_topic_name)
{
  if (ros::param::has("/imu_topic"))
  {
    ros::param::get("/imu_topic", imu_topic_name);
  }
  else
  {
    ROS_WARN("No IMU topic set in yaml file");
  }
}

void setDVLTopicNameFromYawl(std::string& dvl_topic_name)
{
  if (ros::param::has("/dvl_topic"))
  {
    ros::param::get("/dvl_topic", dvl_topic_name);
  }
  else
  {
    ROS_WARN("No DVL topic set in yaml file");
  }
}

void setPressureZTopicNameFromYaml(std::string& pressure_Z_topic_name)
{
  if (ros::param::has("/pressureZ_topic"))
  {
    ros::param::get("/pressureZ_topic", pressure_Z_topic_name);
  }
  else
  {
    ROS_WARN("No PressureZ topic set in yaml file");
  }
}

void setPublishrateFromYaml(int& publish_rate)
{
  if (ros::param::has("/publish_rate"))
  {
    ros::param::get("/publish_rate", publish_rate);
  }
  else
  {
    ROS_WARN("No publish rate set, using default: %i ", publish_rate);
  }
}

void setRdvlFromYamlFile(Matrix3d& R_dvl)
{
  XmlRpc::XmlRpcValue R_dvlConfig;
  if (ros::param::has("/R_dvl"))
  {
    ros::param::get("/R_dvl", R_dvlConfig);
    int matrix_size = R_dvl.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      for (int j = 0; j < matrix_size; j++)
      {
        std::ostringstream ostr;
        ostr << R_dvlConfig[matrix_size * i + j];
        std::istringstream istr(ostr.str());
        istr >> R_dvl(i, j);
      }
    }
  }
  else
  {
    ROS_FATAL("No measurement covariance for DVL (R_dvl) set in parameter file");
    ROS_BREAK();
  }
}

void setRpressureZFromYamlFile(Matrix<double, 1, 1>& R_pressureZ)
{
  XmlRpc::XmlRpcValue R_pressureZConfig;

  if (ros::param::has("/R_pressureZ"))
  {
    ros::param::get("/R_pressureZ", R_pressureZConfig);
    int matrix_size = R_pressureZ.rows();

    for (int i = 0; i < matrix_size; i++)
    {
      std::ostringstream ostr;
      ostr << R_pressureZConfig[i];
      std::istringstream istr(ostr.str());
      istr >> R_pressureZ(i);
    }
  }
  else
  {
    ROS_FATAL("No measurement covariance for pressure sensor (R_pressureZ) set in parameter file!");
    ROS_BREAK();
  }
}


int main(int argc, char* argv[])
{
  ros::init(argc, argv, "eskf");
  ros::NodeHandle nh;
  ros::NodeHandle pnh("~");
  ESKF_Node eskf_node(nh, pnh);
  ros::spin();
  return 0;
}