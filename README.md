## CCO-VOXEL: Chance Constrained Optimization over Uncertain Voxel-Grid Representation for Safe Trajectory Planning

CCO VOXEL is an algorithm that computes trajectory plans with probabilistic safety guarantees in real-time directly on the voxel-grid representation of the world. CCO-VOXEL maps the distribution over the distance to the closest obstacle to a distribution over collision-constraint violation and computes an optimal trajectory that minimizes the violation probability.

![](https://github.com/sudarshan-s-harithas/CCO-VOXEL/blob/main/Images/teaser.png?raw=true)
 

#### Preprint: https://arxiv.org/abs/2110.02904 

## Setup 

Our evaluations were done on an 8-Core Intel Core i7-10870H processor, with 16GB RAM and 500GB HDD, running Ubunut 20.04 and ros noetic. we do recommend using powerfull setup. 

### Pre-requisites

[ROS Noetic](http://wiki.ros.org/noetic/Installation/Ubuntu) <br />
[Octomap](http://wiki.ros.org/octomap) <br />
[Mavros](https://docs.px4.io/master/en/ros/mavros_installation.html) <br />
[PX4-Autopilot](https://docs.px4.io/master/en/dev_setup/dev_env_linux_ubuntu.html#gazebo-jmavsim-and-nuttx-pixhawk-targets)<br />
[QGroundControl](https://docs.qgroundcontrol.com/master/en/getting_started/download_and_install.html) (Optional)<br />

### Build
```
cd ~/catkin_ws/src/ 
git clone https://github.com/sudarshan-s-harithas/CCO-VOXEL.git
catkin_make -j4
source devel/setup.bash 
```

### Evaluation 

We provide Gazebo environments to test our algorithm that can be found [here](https://github.com/sudarshan-s-harithas/CCO-VOXEL/tree/main/CCO_VOXEL/worlds), Please follow the the instructions given [here](https://github.com/sudarshan-s-harithas/CCO-VOXEL/tree/main/CCO_VOXEL#origanization-of-your-working-directories) before continuing with the execution of the programs.    

![](https://github.com/sudarshan-s-harithas/CCO-VOXEL/blob/main/Images/simulation.gif)

Run the following commands to start CCO-VOXEL <br />

First within the /PX4-Autopilot folder run the commands below to start SITL 

```
Terminal1: 
source Tools/setup_gazebo.bash $(pwd) $(pwd)/build/px4_sitl_default
export ROS_PACKAGE_PATH=$ROS_PACKAGE_PATH:$(pwd)
export ROS_PACKAGE_PATH=$ROS_PACKAGE_PATH:$(pwd)/Tools/sitl_gazebo
roslaunch px4 mavros_posix_sitl.launch
```
Once SITL is started run Mapping_noise.launch to start *OctoMap* and Rviz, the *noise* parameter can be changed in the launch file to update the added noise levels into the point cloud.  
```
Terminal2: 
roslaunch CCO_VOXEL Mapping_noise.launch
```
After starting the mapping process, takeoff the drone either using QGround Control or through the *commander takeoff* command, and start the planner and controller,  and update the parameter *path_to_weights* to point the weights.csv file in the *CCO_VOXEL_Planner.launch* file. 
```
Terminal3: 
roslaunch CCO_VOXEL CCO_VOXEL_Planner.launch

Terminal4: 

rosrun CCO_VOXEL Controller
```

Once the programs are running use the *2D Nav Goal* tool of rviz to set the goal point. 

### PlayGround

For better understanding of the algorithm please check the codes here[https://github.com/sudarshan-s-harithas/CCO-VOXEL/tree/main/CCO_VOXEL/PlayGround]

### Acknowledgements 
Our code is built upon [Fast-Planner](https://github.com/HKUST-Aerial-Robotics/Fast-Planner), we use their front end implementation of the *kinodynamic A* * with a difference that we use MMD as a part of the Edge Cost that connects two nodes of the graph. 
