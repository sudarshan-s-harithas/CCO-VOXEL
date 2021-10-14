#include"utils.h"


namespace Visualizer
{


	class Trajectory_visualizer
	{

	public:

		int visulize_sampled_trajectories(Eigen::MatrixXd TopX, Eigen::MatrixXd TopY, Eigen::MatrixXd TopZ , int num_samples , int num_pts_per_traj , ros::Publisher sample_trajectory_pub  );



	};
} 


int Visualizer::Trajectory_visualizer::visulize_sampled_trajectories(Eigen::MatrixXd TopX, Eigen::MatrixXd TopY, Eigen::MatrixXd TopZ, int num_samples , int num_pts_per_traj, ros::Publisher sample_trajectory_pub  )
{

int color_counter = 0;
for( int i =0 ; i < num_samples ; i++){
color_counter++;
visualization_msgs::Marker marker;
marker.type = visualization_msgs::Marker::LINE_LIST;
marker.action = visualization_msgs::Marker::ADD;

geometry_msgs::Point pt;

for( int j =0 ; j< num_pts_per_traj ; j++){

// visualization_msgs::Marker marker;
marker.header.frame_id = "map";
marker.header.stamp = ros::Time();
marker.ns = "my_namespace";
marker.id = color_counter;


pt.x = TopX( i, j );
pt.y  = TopY( i, j );
pt.z =TopZ( i, j );

marker.points.push_back(pt);


marker.scale.x = 0.1;
marker.scale.y = 0.1;
marker.scale.z = 0.1;
marker.color.a = 1; 
marker.color.r = 1.0 ; //(num_samples - float(color_counter))/num_samples;
marker.color.g = 0.5;
marker.color.b = 0.2*(num_samples - float(color_counter))/num_samples;

// std::cout << TopX( i, j ) << "  " << TopY( i, j ) << "  " << TopZ( i, j ) << std::endl;

}

sample_trajectory_pub.publish( marker);
} 


return 0;

}