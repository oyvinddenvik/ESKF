#pragma once

#include "ESKF.h"
#include "common.h"

#include "ros/ros.h"
#include "sensor_msgs/Imu.h"
#include "nav_msgs/Odometry.h"


class ESKF_Node
{
public:
    ESKF_Node();
    ~ESKF_Node();

    void imuCallback(const sensor_msgs::Imu::ConstPtr& imu_Message_data);


private:

    ESKF eskf;
    
};




