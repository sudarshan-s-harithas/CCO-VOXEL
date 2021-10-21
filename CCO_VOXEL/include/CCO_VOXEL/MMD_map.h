#include<dynamicEDT3D/dynamicEDTOctomap.h> 



namespace MMD_Map{

	class MMD_Map_Functions{

	public:

		DynamicEDTOctomap* EDT_Map; 
		Eigen::Matrix< float , 100, 5> Weights;
		octomap::OcTree *tree = new octomap::OcTree(0.05);
		octomap::point3d start, end;
		double octree_resolution = 0.2; 
		double resolution_factor = 1/octree_resolution ;
		octomap::OcTree *octree_for_mmd;
		octomap::point3d mmd_map_start , mmd_map_end;

		int key_x , key_y , key_z ; 

		// void set_edt_map( DynamicEDTOctomap *ptr);

		float MMD_cost_per_point_temp , mean_obstacle_distance_temp; 
		int xkey_temp , ykey_temp , zkey_temp;

        struct MMD_Cost{
            float MMD_cost_per_point;
            float mean_obstacle_distance;
            // octomap::point3d  mean_obstacle_location; 
        };

        struct MMD_Map_key{

            int xkey;
            int ykey;
            int zkey;
            bool operator<(const MMD_Map_key& t) const
            {
                return (this->xkey < t.xkey);
            }
        };


        struct Occupancy_Map_key{

            int key_valx;
            int key_valy;
            int key_valz;
            bool operator<(const Occupancy_Map_key& t) const
            {
                return (this->key_valx < t.key_valx);
            }
        };



        std::map<MMD_Map_key , MMD_Cost> MMD_data;
        std::map<Occupancy_Map_key , bool> Occupancy_data;

        // void update_mmd_map();

        void update_MMD_Map(  DynamicEDTOctomap* ptr , visualization_msgs::MarkerArray mdd_marker ,  ros::Publisher MMD_map_pub);
        double compute_MMD_linear_transforms(Eigen::MatrixXf actual_distribution);
        void assign_weights_for_MMD();
        void convert_point_to_key(octomap::point3d inpt  ,  int &key_x , int &key_y , int &key_z  );
        std::tuple< int , int , int> convert_point_to_key_external( octomap::point3d inPt);
        float get_MMD_cost_per_point( octomap::point3d Query_Point);
        double  compute_EDT_interpolation( float distance_at_query_point);

        std::tuple< bool, float , bool, float , bool> Point_within_map_and_mmd_value( octomap::point3d  state_pos_start ,octomap::point3d  state_pos_end );
        bool get_occupancy(octomap::OcTreeKey Node_key );


    private:

    	void convert_point_to_key_tuple();

    	octomap::point3d convertToOctomapPt(Eigen::Vector3d pt);
    	Eigen::Vector3d convertToEigenPt(octomap::point3d pt); 





	};


}



// void MMD_Map::MMD_Map_Functions::update_mmd_map()
// {

// MMD_Cost MMD_c;

// MMD_c.MMD_cost_per_point = MMD_cost_per_point_temp;
// MMD_c.mean_obstacle_distance = mean_obstacle_distance_temp;


// MMD_Map_key MMK_new;


// MMK_new.xkey = xkey_temp;
// MMK_new.ykey = ykey_temp;
// MMK_new.zkey = zkey_temp;


 
// MMD_Map[MMK_new] = MMD_c;
// }

double MMD_Map::MMD_Map_Functions::compute_EDT_interpolation( float distance_at_query_point )
{

    Eigen::MatrixXf coefficient_matrix(1,7);
    Eigen::MatrixXf distance_matrix(1,7);
    Eigen::MatrixXf result(1 ,1);

    coefficient_matrix(0,0) = 7720.61614034;
    coefficient_matrix(0,1) =-12367.93260857;
    coefficient_matrix(0,2) =-1541.47268526;
    coefficient_matrix(0,3) = 14365.74417541 ;
    coefficient_matrix(0,4) =-10918.17820688;
    coefficient_matrix(0,5) = 3339.69565025;
    coefficient_matrix(0,6) = -373.83984028;

    for( int i=0 ; i< 7 ; i++){
        distance_matrix(0 ,i) = pow( distance_at_query_point ,i);
    }

    result = coefficient_matrix*distance_matrix.transpose();

    return result(0,0);

}



std::tuple< bool, float , bool, float , bool>  MMD_Map::MMD_Map_Functions::Point_within_map_and_mmd_value( octomap::point3d  state_pos_start ,octomap::point3d  state_pos_end ){


std::tuple< int , int , int> map_keys_start ; 
std::tuple< int , int , int>  map_keys_end ;

map_keys_start = convert_point_to_key_external( state_pos_start);
map_keys_end = convert_point_to_key_external( state_pos_end);
bool start_point_in_map = true;
bool end_point_in_map = true;


std::tuple< bool, float , bool, float , bool > point_and_mmd;



MMD_Map_key Query_key_data_start;

Query_key_data_start.xkey = std::get<0>(map_keys_start);
Query_key_data_start.ykey = std::get<1>(map_keys_start);
Query_key_data_start.zkey = std::get<2>(map_keys_start);

if( MMD_data.find(Query_key_data_start) == MMD_data.end()){
    std::get<0>(point_and_mmd)  = false; 
    start_point_in_map = false; 
}
else{

    std::get<0>(point_and_mmd)  = true; 
    start_point_in_map = true; 
}


MMD_Map_key Query_key_data_end;


Query_key_data_end.xkey = std::get<0>(map_keys_end);
Query_key_data_end.ykey = std::get<1>(map_keys_end);
Query_key_data_end.zkey = std::get<2>(map_keys_end);


if( MMD_data.find(Query_key_data_end) == MMD_data.end()){
    std::get<2>(point_and_mmd)  = false; 
    end_point_in_map = false;
}
else{
    std::get<2>(point_and_mmd) = true; 
    end_point_in_map = true;
}



if( start_point_in_map){

float MMD_value_start;

MMD_value_start = MMD_data[Query_key_data_start].MMD_cost_per_point ;

 std::get<1>(point_and_mmd) = MMD_value_start;

}


if( end_point_in_map){

float MMD_value_end;

MMD_value_end = MMD_data[Query_key_data_end].MMD_cost_per_point ;
 std::get<3>(point_and_mmd) = MMD_value_end;

}

bool to_continue = true;

if( end_point_in_map && start_point_in_map){

    to_continue = true;
}
else{
    to_continue = false;
}

std::get<4>(point_and_mmd) = to_continue;

return point_and_mmd; 

}


void  MMD_Map::MMD_Map_Functions::convert_point_to_key( octomap::point3d  inpt , int &key_x , int &key_y , int &key_z )
{

    Eigen::Vector3d temp_key_val  ;
    temp_key_val = convertToEigenPt( inpt);

    double xin , yin ,zin;

    xin = double( temp_key_val.x() );
    yin = double ( temp_key_val.y());
    zin = double( temp_key_val.z());

    key_x = (int)floor(xin*resolution_factor);
    key_y = (int)floor(yin*resolution_factor);
   key_z = (int)floor(zin*resolution_factor);



}



std::tuple< int , int , int> MMD_Map::MMD_Map_Functions::convert_point_to_key_external( octomap::point3d inPt) 
{

    Eigen::Vector3d temp_key_val  ;
    temp_key_val = convertToEigenPt( inPt);

    std::tuple< int, int , int > tuple_key_values;

    double xin , yin ,zin;

    xin = double( temp_key_val.x() );
    yin = double ( temp_key_val.y());
    zin = double( temp_key_val.z());

    std::get<0>(tuple_key_values) = (int)floor(xin*resolution_factor);
    std::get<1>(tuple_key_values) = (int)floor(yin*resolution_factor);
   std::get<2>(tuple_key_values) = (int)floor(zin*resolution_factor);

return tuple_key_values;

}



void MMD_Map::MMD_Map_Functions::assign_weights_for_MMD(){


	int rows = 100;
int cols =5;
 std::string file ="../src/CCO_VOXEL/include/CCO_VOXEL/weight.csv" ;
  std::ifstream in(file);
  
  std::string line;

  int row = 0;
  int col = 0;

  Eigen::MatrixXf res = Eigen::MatrixXf(rows, cols);

  if (in.is_open()) {

    while (std::getline(in, line)) {

      char *ptr = (char *) line.c_str();
      int len = line.length();

      col = 0;

      char *start = ptr;
      for (int i = 0; i < len; i++) {

        if (ptr[i] == ',') {
          res(row, col++) = atof(start);
          start = ptr + i + 1;
        }
      }
      res(row, col) = atof(start);

      row++;
    }

    in.close();
  }

    Weights = res;

  std::cout << Weights.rows() << " " << Weights.cols() << std::endl;
}



double MMD_Map::MMD_Map_Functions::compute_MMD_linear_transforms(Eigen::MatrixXf actual_distribution)
{

	Eigen::MatrixXf transformed_features(1,  5);


transformed_features = actual_distribution*Weights;

// std::cout << transformed_features << std::endl;

int num_samples_of_distance_distribution =transformed_features.cols() ;
Eigen::MatrixXf squared_dist(1, num_samples_of_distance_distribution);
Eigen::MatrixXf temp(1, num_samples_of_distance_distribution);
Eigen::MatrixXf alpha_weights(1, num_samples_of_distance_distribution);

alpha_weights.setOnes();


// temp = actual_distribution.array();

squared_dist = (transformed_features.array()).square();

Eigen::MatrixXf kernel_matrix1(num_samples_of_distance_distribution, num_samples_of_distance_distribution);
Eigen::MatrixXf temp_kernel_matrix1(num_samples_of_distance_distribution, num_samples_of_distance_distribution);
Eigen::MatrixXf temp_kernel_matrix2(num_samples_of_distance_distribution, num_samples_of_distance_distribution);


Eigen::MatrixXf One_matrix( num_samples_of_distance_distribution, num_samples_of_distance_distribution);
One_matrix.setOnes();


// std::cout << " Entered the main funciton /////////////" << std::endl;
temp_kernel_matrix1 = transformed_features.transpose()*transformed_features;

temp_kernel_matrix2 = squared_dist.transpose()*squared_dist;

kernel_matrix1 = One_matrix + 2*temp_kernel_matrix1 + temp_kernel_matrix2;

Eigen::Matrix<float, 1,1>  res; 
Eigen::Matrix<float, 1,1>  res1; 
Eigen::Matrix<float, 1,1>  res2; 
Eigen::Matrix<float, 1,1>  res3; 

float temp_value =0; 

     

double res_val;
res1 = alpha_weights*kernel_matrix1*(alpha_weights.transpose()) ; 

res2(0 , 0 ) = num_samples_of_distance_distribution;
res3(0,0) = num_samples_of_distance_distribution ;

res(0,0) = res1(0,0) - 2*res2(0,0) +res3(0,0) ;
res_val = double( res(0,0) );

return res_val; 


}



float MMD_Map::MMD_Map_Functions::get_MMD_cost_per_point( octomap::point3d Query_Point){

    std::tuple< int, int, int > map_keys; 

 // auto start =  std::chrono::high_resolution_clock::now();

map_keys = convert_point_to_key_external( Query_Point);

MMD_Map_key Query_key_data;

Query_key_data.xkey = std::get<0>(map_keys);
Query_key_data.ykey = std::get<1>(map_keys);
Query_key_data.zkey = std::get<2>(map_keys);

float MMD_value;

MMD_value = MMD_data[Query_key_data].MMD_cost_per_point ;

// auto stop =  std::chrono::high_resolution_clock::now();

// auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);



// std::cout << " --------------------------- " <<    MMD_value << "   "  << duration.count()  << std::endl; 

return MMD_value;


}


void MMD_Map::MMD_Map_Functions::update_MMD_Map(  DynamicEDTOctomap* ptr , visualization_msgs::MarkerArray mdd_marker ,  ros::Publisher MMD_map_pub)
{

   unsigned char maxDepth = 16;
    uint32_t shape = visualization_msgs::Marker::CUBE;
    int count  = 0;
    octomap::OcTreeKey key_from_octomap;
    int num_samples_of_distance_distribution =100;
     std::random_device rd{};
     std::mt19937 gen{rd()};

     std::normal_distribution<> noise{0,2};    
     Eigen::MatrixXf noise_distribution(1, num_samples_of_distance_distribution);
     Eigen::MatrixXf radius(1, num_samples_of_distance_distribution);  
     // radius.setOnes();
     radius = Eigen::MatrixXf::Constant( 1, num_samples_of_distance_distribution, 0.75); 
     Eigen::MatrixXf actual_distance( 1, num_samples_of_distance_distribution);
     Eigen::MatrixXf actual_distribution( 1, num_samples_of_distance_distribution);
     Eigen::MatrixXf zero_matrix( 1, num_samples_of_distance_distribution);
     zero_matrix.setZero();
     assign_weights_for_MMD();
     int mmd_num_counter = 0 ;


    for (int i=0 ; i<num_samples_of_distance_distribution; i++){


      noise_distribution(0 ,i) = radius(0 ,i) - float(noise(gen));
    }



    for(octomap::OcTree::leaf_bbx_iterator it = octree_for_mmd->begin_leafs_bbx(mmd_map_start, mmd_map_end, maxDepth), bbx_end = octree_for_mmd->end_leafs_bbx(); it != bbx_end; std::advance(it,3))
    {




        // std::cout<<it.getCoordinate()<<std::endl;
        octomap::point3d pt = it.getCoordinate();
        // Eigen::Matrix<int , 3,1,> temp_key_vector;

        convert_point_to_key( pt  , key_x , key_y , key_z);

        MMD_Map_key Key_data;

        Key_data.xkey= key_x;
        Key_data.ykey =key_y;
        Key_data.zkey=key_z;

        Occupancy_Map_key OCC_keys;



        key_from_octomap = it.getKey();  // conversion from MMD_Map_key to key_from_octomap is key_from_octomap = MMD_Map_key + max_val_tree the mac_val_tree is a constant

        bool occ_info = octree_for_mmd->isNodeOccupied(*it);

        


        OCC_keys.key_valx =  key_from_octomap[0];
        OCC_keys.key_valy =  key_from_octomap[1];
        OCC_keys.key_valz =  key_from_octomap[2];

        Occupancy_data[OCC_keys] = occ_info;

        if(!ros::ok())
        {
            break;
        }

        visualization_msgs::Marker marker;
        marker.header.frame_id = "map";
        marker.header.stamp = ros::Time::now();
        marker.ns = count;//pt.x() + pt.y() + pt.z();
        marker.id = count;
        marker.type = shape;

        count++;

        // set the marker action. Options are ADD and DELETE
        marker.action = visualization_msgs::Marker::ADD;

        // set the pose of the marker with respect to the global frame of reference
        marker.pose.position.x = pt.x();
        marker.pose.position.y = pt.y();
        marker.pose.position.z = pt.z();

        marker.pose.orientation.x = 0;
        marker.pose.orientation.y = 0;
        marker.pose.orientation.z = 0;
        marker.pose.orientation.w = 1;

        // set the scale of the marker
        marker.scale.x = 0.30;
        marker.scale.y = 0.30;
        marker.scale.z = 0.30;

        // set the color of the cell based on the distance from the obstacle
        float dist = ptr->getDistance(pt);
        float MMD_val =0;

        double threshold_val = 0.75 - ( dist -2);

        if( threshold_val > 0){
            actual_distance = Eigen::MatrixXf::Constant( 1, num_samples_of_distance_distribution, dist);
            actual_distribution = zero_matrix.cwiseMax( noise_distribution - actual_distance );

            MMD_val = compute_MMD_linear_transforms(actual_distribution ); // compute_MMD_linear_transforms(actual_distribution ); //get_EDT_cost(actual_distribution );  //compute_EDT_interpolation(dist ); //get_EDT_cost(dist ); *

            mmd_num_counter++;
            // std::cout <<  mmd_num_counter << std::endl;
    }

    if(threshold_val <= 0){

         

        MMD_val =0;
    }

    unsigned short int a;

            //std::cout << MMD_val <<  " MMD_val " << std::endl;

           MMD_Cost MMDc; 

            MMDc.MMD_cost_per_point = MMD_val ;
            MMDc.mean_obstacle_distance = dist;

             MMD_data[Key_data]= MMDc;

             // std::cout << MMD_data[Key_data].MMD_cost_per_point <<  " " << MMD_val << std::endl;

             // std::cout << Key_data.xkey << " " << Key_data.ykey << "  " << Key_data.zkey << std::endl;

             // std::cout << key_from_octomap[0] << "  "  <<  key_from_octomap[1] << " "<<  key_from_octomap[2] <<  std::endl;
              // std::cout << "--------------" << std::endl;


        marker.color.r =  std::sqrt(std::sqrt(MMD_val/10000.0));
        marker.color.g = (1 - MMD_val/10000.0)*(1 - MMD_val/10000.0);
        marker.color.b = 0.0;
        marker.color.a = (1 - dist/1000.0)*(1 - dist/1000.0);
                    
        marker.lifetime = ros::Duration();
        mdd_marker.markers.push_back(marker);



    }

    MMD_map_pub.publish(mdd_marker);
    // delta_Vector.x() = gradient_stepx;
    // delta_Vector.y() = gradient_stepy;
    // delta_Vector.z() = gradient_stepz;
    // delta = convertToOctomapPt(delta_Vector);
    // Eigen::Vector3d gradient_Vector_per_point;
/*
for(octomap::OcTree::leaf_bbx_iterator grad_it = tree->begin_leafs_bbx(start, end, maxDepth), bbx_end = tree->end_leafs_bbx(); grad_it != bbx_end; std::advance(grad_it,3))
    {

        MMD_Gradients upload_grad;


        octomap::point3d grad_pt = grad_it.getCoordinate();

        // std::cout << grad_pt.x() << std::endl;
        convert_point_to_key( grad_pt  , key_x , key_y , key_z);
        gradient_Vector_per_point = compute_gradients_MMD_map( grad_pt );

        upload_grad.gradient_x = float( gradient_Vector_per_point.x());
        upload_grad.gradient_y = float( gradient_Vector_per_point.y());
        upload_grad.gradient_z = float( gradient_Vector_per_point.z());

        MMD_grad_key upload_Key_data;

        upload_Key_data.x_grad_key = key_x;
        upload_Key_data.y_grad_key = key_y;
        upload_Key_data.z_grad_key = key_z;

        MMD_gradient_data[upload_Key_data]= upload_grad;

// std::cout << "************" << std::endl;


//  std::cout << "  " << gradient_Vector_per_point.x() << "  " << gradient_Vector_per_point.y() << " " << gradient_Vector_per_point.z() << std::endl;

// std::cout << "--------------" << std::endl;

    }

*/


}



octomap::point3d MMD_Map::MMD_Map_Functions::convertToOctomapPt(Eigen::Vector3d pt)
{
    octomap::point3d p(pt(0), pt(1), pt(2));
    return p;
}


////////////////////////////////////////////////////////////////////////////////////////
Eigen::Vector3d MMD_Map::MMD_Map_Functions::convertToEigenPt(octomap::point3d pt)
{
    Eigen::Vector3d p(pt.x(), pt.y(), pt.z());
    return p;
}




bool MMD_Map::MMD_Map_Functions::get_occupancy(octomap::OcTreeKey Node_key )
{
Occupancy_Map_key OMK;

OMK.key_valx = Node_key[0];
OMK.key_valy = Node_key[1];
OMK.key_valz = Node_key[2];

bool occ_res =false ;
occ_res = Occupancy_data[OMK];

return occ_res;
}