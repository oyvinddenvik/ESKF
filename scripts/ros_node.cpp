#include"common.h"
#include"ESKF.h"
#include<Eigen/Core>
#include"ros_node.h"


/*
void ESKF_Node::getParametersFromYamlFile()
{




}
*/



ESKF_Node::ESKF_Node(const ros::NodeHandle& nh, const ros::NodeHandle& pnh) : nh_{pnh}, init_{false}, eskf_{R_ACC,R_ACCBIAS,R_GYRO,R_GYROBIAS,P_GYRO_BIAS,P_ACC_BIAS,S_A,S_G,S_DVL,S_INC}
{
    int publish_rate{125}; // Change this to include params
    ROS_INFO("Subscribing to IMU /manta/imu");
    // Subscribe to IMU
    subscribeIMU_ = nh_.subscribe<sensor_msgs::Imu>("/manta/imu",1000,&ESKF_Node::imuCallback,this, ros::TransportHints().tcpNoDelay(true));


    ROS_INFO("Publishing State");
    publishPose_ = nh_.advertise<nav_msgs::Odometry>("pose",1);
    
    
    pubTImer_= nh_.createTimer(ros::Duration(1.0f/publish_rate), &ESKF_Node::publishPoseState, this);
}








// IMU Subscriber
void ESKF_Node::imuCallback(const sensor_msgs::Imu::ConstPtr& imu_Message_data)
{
    double Ts{0};
    int imu_publish_rate{DEFAULT_IMU_RATE};
    Vector3d raw_acceleration_measurements = Vector3d::Zero();
    Vector3d raw_gyro_measurements = Vector3d::Zero();

    Ts = (1.0/imu_publish_rate);
    raw_acceleration_measurements << imu_Message_data->linear_acceleration.x,
                                   imu_Message_data->linear_acceleration.y,
                                   imu_Message_data->linear_acceleration.z;

    raw_gyro_measurements << imu_Message_data->angular_velocity.x,
                             imu_Message_data->angular_velocity.y,
                             imu_Message_data->angular_velocity.z;

    ROS_INFO("Acceleration_x: %f",imu_Message_data->linear_acceleration.x);
    eskf_.predict(raw_acceleration_measurements,raw_gyro_measurements,Ts);       
              
}


void ESKF_Node::publishPoseState(const ros::TimerEvent&)
{
    nav_msgs::Odometry odom_msg;
    static size_t trace_id{0};
    const VectorXd pose = eskf_.getPose();

    odom_msg.header.frame_id = "/eskf_link";
    odom_msg.header.seq = trace_id++;
    odom_msg.header.stamp = ros::Time::now();
    odom_msg.pose.pose.position.x = pose(StateMemberX);
    odom_msg.pose.pose.position.y = pose(StateMemberY);
    odom_msg.pose.pose.position.z = pose(StateMemberZ);
    odom_msg.twist.twist.linear.x = pose(StateMemberVx);
    odom_msg.twist.twist.linear.y = pose(StateMemberVy);
    odom_msg.twist.twist.linear.z = pose(StateMemberVz);
    odom_msg.pose.pose.orientation.w = pose(StateMemberQw);
    odom_msg.pose.pose.orientation.x = pose(StateMemberQx);
    odom_msg.pose.pose.orientation.y = pose(StateMemberQy);
    odom_msg.pose.pose.orientation.z = pose(StateMemberQz);

    //ROS_INFO("StateX: %f",odom_msg.pose.pose.position.x);
    publishPose_.publish(odom_msg);

}

