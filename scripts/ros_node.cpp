#include"common.h"
#include"ESKF.h"
#include<Eigen/Core>
#include"ros_node.h"







// IMU Subscriber
void ESKF_Node::imuCallback(const sensor_msgs::Imu::ConstPtr& imu_Message_data)
{

    Vector3d raw_acceration_measurements = Vector3d::Zero();
    Vector3d raw_gyro_measurements = Vector3d::Zero();

    raw_acceration_measurements << imu_Message_data->linear_acceleration.x,
                                   imu_Message_data->linear_acceleration.y,
                                   imu_Message_data->linear_acceleration.z;
    
    raw_gyro_measurements << imu_Message_data->angular_velocity.x,
                             imu_Message_data->angular_velocity.y,
                             imu_Message_data->angular_velocity.z;

}

