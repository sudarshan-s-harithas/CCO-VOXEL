#include"utils.h"
#include <random>
#include<math.h>
#include<dynamicEDT3D/dynamicEDTOctomap.h>
#include"Bernstein_test.h"

Bernstein::Bernstein_polynomial_function BPF; 


namespace Optimization
{

	class CEOptim{


	public:

		std::vector<Eigen::Vector3d> Bernstein_coeff;
		std::vector<Eigen::Vector3d> perturbed_Bernstein_coeff;


		nav_msgs::Path optimize_generated_path(std::vector<Eigen::Vector3d>  currTraj , DynamicEDTOctomap *edt_map_ptr , float time_to_desination   );


	};






}



nav_msgs::Path Optimization::CEOptim::optimize_generated_path(std::vector<Eigen::Vector3d>  currTraj , DynamicEDTOctomap *edt_map_ptr , float time_to_desination  )
{

float distance;
octomap::point3d pt; 
pt.x() = 5.0;
pt.y() = 4.0;
pt.z() = 3.5;

distance = edt_map_ptr->getDistance(pt);
int num_points_path =0 ;

for(auto i = currTraj.begin(); i != currTraj.end(); i++)
{

num_points_path = num_points_path+1;

}
nav_msgs::Path refined_path;

Eigen::MatrixXd Path_vector( num_points_path , 3);
int count=0;
for(auto i = currTraj.begin(); i != currTraj.end(); i++)
{

Eigen::Vector3d pos = *i;
Path_vector(count ,0  ) = pos(0);
Path_vector(count ,1  ) = pos(1);
Path_vector(count ,2  ) = pos(2);
count++;


}


std::cout<< " done this" <<std::endl;

Bernstein_coeff  = BPF.get_coefficients( currTraj , num_points_path ,time_to_desination );
perturbed_Bernstein_coeff = BPF.perturb_coefficients(Bernstein_coeff);
return refined_path;


}