## CCO-VOXEL: Chance Constrained Optimization over Uncertain Voxel-Grid Representation for Safe Trajectory Planning

CCO VOXEL is an algorithm that computes trajectory plans with probabilistic safety guarantees in real-time directly on the voxel-grid representation of the world. CCO-VOXEL maps the distribution over the distance to the closest obstacle to a distribution over collision-constraint violation and computes an optimal trajectory that minimizes the violation probability.
 
This page is under active development  


## Setup 

Our evaluations were done on an 8-Core Intel Core i7-10870H processor, with 16GB RAM and 500GB HDD, running Ubunut 20.04 and ros noetic. we do recommend using powerfull setup. 

### Pre-requisites

[ROS Noetic](http://wiki.ros.org/noetic/Installation/Ubuntu) <br />
[Octomap](http://wiki.ros.org/octomap) <br />
[Mavros](https://docs.px4.io/master/en/ros/mavros_installation.html) <br />
[QGroundControl](https://docs.qgroundcontrol.com/master/en/getting_started/download_and_install.html) (Optional)<br />
