#pragma once

#include "ESKF.h"
#include "common.h"

#include "ros/ros.h"
#include "sensor_msgs/Imu.h"
#include "nav_msgs/Odometry.h"





const Eigen::Matrix3d R_ACC((Eigen::Matrix3d() << 4,0,0,0,4,0,0,0,4).finished());
const Eigen::Matrix3d R_ACCBIAS((Eigen::Matrix3d() << 6e-5,0,0,0,6e-5,0,0,0,6e-5).finished());
const Eigen::Matrix3d R_GYRO((Eigen::Matrix3d() << 12e-3,0,0,0,12e-3,0,0,0,12e-3).finished());
const Eigen::Matrix3d R_GYROBIAS((Eigen::Matrix3d() << 3e-7,0,0,0,3e-7,0,0,0,3e-7).finished());
constexpr double P_GYRO_BIAS{0.0001};
constexpr double P_ACC_BIAS{0.0001};
const Eigen::Matrix3d S_A((Eigen::Matrix3d() << -0.9990,1.9804e-04,-0.0450,1.3553e-20,1.0,0.0044,0.0450,0.0044,-0.9990).finished());
const Eigen::Matrix3d S_G((Eigen::Matrix3d() << -0.9990,1.9804e-04,-0.0450,1.3553e-20,1.0,0.0044,0.0450,0.0044,-0.9990).finished());
const Eigen::Matrix3d S_DVL((Eigen::Matrix3d() << 0.0001,0,0,0,0.0001,0,0,0,0.0001).finished());
const Eigen::Matrix3d S_INC((Eigen::Matrix3d() << 0.0001,0,0,0,0.0001,0,0,0,0.0001).finished());
//const Eigen::VectorXd INITIAL_NOMINAL_STATE((Eigen::VectorXd() << 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15).finished());
//const Eigen::VectorXd INITIAL_NOMINAL_STATE = (Eigen::VectorXd(16) << 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16).finished();
//const Eigen::MatrixXd Initial_P =(Eigen::MatrixXd(2,2) << 1,2,3,4).finished();






class ESKF_Node
{
public:
    ESKF_Node(const ros::NodeHandle& nh, const ros::NodeHandle& pnh);
    //~ESKF_Node();

    
    //void getParametersFromYamlFile();


private:
    ros::NodeHandle nh_;
    bool init_;

    
    //VectorXd initialNominalState_;
    //Matrix<double,ERROR_STATE_SIZE,ERROR_STATE_SIZE> initialP_;
   
    ESKF eskf_;



    // Load from Yaml file
    parametersInESKF loadParametersFromYamlFile();

    // ROS subscribers
    ros::Subscriber subscribeIMU_;
    ros::Timer pubTImer_;


    // ROS publisher
    ros::Publisher publishPose_;
    

    // Callbacks
    void imuCallback(const sensor_msgs::Imu::ConstPtr& imu_Message_data);
    void publishPoseState(const ros::TimerEvent&);

};




