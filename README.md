## CCO-VOXEL: Chance Constrained Optimization over Uncertain Voxel-Grid Representation for Safe Trajectory Planning

CCO VOXEL is an algorithm that computes trajectory plans with probabilistic safety guarantees in real-time directly on the voxel-grid representation of the world. CCO-VOXEL maps the distribution over the distance to the closest obstacle to a distribution over collision-constraint violation and computes an optimal trajectory that minimizes the violation probability.

![](https://github.com/sudarshan-s-harithas/CCO-VOXEL/blob/main/Images/teaser.png?raw=true)
 
#### This page is under active development  

#### Preprint: https://arxiv.org/abs/2110.02904 

## Setup 

Our evaluations were done on an 8-Core Intel Core i7-10870H processor, with 16GB RAM and 500GB HDD, running Ubunut 20.04 and ros noetic. we do recommend using powerfull setup. 

### Pre-requisites

[ROS Noetic](http://wiki.ros.org/noetic/Installation/Ubuntu) <br />
[Octomap](http://wiki.ros.org/octomap) <br />
[Mavros](https://docs.px4.io/master/en/ros/mavros_installation.html) <br />
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


```
Terminal1: 
source Tools/setup_gazebo.bash $(pwd) $(pwd)/build/px4_sitl_default
export ROS_PACKAGE_PATH=$ROS_PACKAGE_PATH:$(pwd)
export ROS_PACKAGE_PATH=$ROS_PACKAGE_PATH:$(pwd)/Tools/sitl_gazebo
roslaunch px4 mavros_posix_sitl.launch

Terminal2: 
roslaunch CCO_VOXEL Mapping_noise.launch

Terminal3: 
rosrun CCO_VOXEL Planner

Terminal4: 

rosrun CCO_VOXEL Controller
```

### Acknowledgements 
TODO
