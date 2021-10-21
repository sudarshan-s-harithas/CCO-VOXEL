/** 
 * Map header file
 * Store the parameters of the map
**/
#pragma once

#include"utils.h"
#include <random>
#include<algorithm>
#include <limits.h>
#include <queue>
#include<tuple> 
#include<math.h>

// MMD_Map::MMD_Map_Functions MMD_costmap;

namespace Map3D
{
    class OctoMapEDT
    {
        public:

            //// data members ////
            octomap::AbstractOcTree *new_tree = NULL; // reads the octomap when available in callback
            double octree_resolution = 0.2; 
            double resolution_factor = 1/octree_resolution ;
            octomap::OcTree *tree = new octomap::OcTree(octree_resolution); // convert to OcTree for EDT calculation            
            octomap::OcTree *MMD_tree = new octomap::OcTree(octree_resolution);    

            octomap::point3d min, max;  // min, max of the map
            octomap::point3d start, end; // start, end of the fixed window EDT
            octomap::point3d mapStart; // location where mapping starts or reset
            float mapUpdateThreshold;  // minimum distance the drone should travel w.r.to previous mapUpdate point for map to update again
            DynamicEDTOctomap* costMap;    // pointer to the EDT of the octomap
            float edtMapRange; // range in x,y,z dimensions (X->forward, y->sideways, z->vertical for the Euclidean Distance Map)
            bool isOctomapUpdated;       // set true when map is updated

            int key_x , key_y , key_z ; 
            double offsetx = 0.08 ;
            double offsety =0.08;
            double offsetz = 0.08 ;

            double gradient_stepx = 1*octree_resolution +offsetx;
            double gradient_stepy = 1*octree_resolution+ offsety;
            double gradient_stepz = 1*octree_resolution+ offsetz;
            Eigen::Vector3d delta_Vector;
            Eigen::Matrix< float , 100, 5> Weights;


            octomap::point3d delta; 
            


            //// member functions /////
            OctoMapEDT();
            bool ifUpdateMap(Eigen::Vector3d pt1);
            void setStartPosition(Eigen::Vector3d pt);
            void setMapRange(Eigen::Vector3d pt);
            bool isInMap(octomap::point3d pt);
            void setMinMax();
            void erase();
            void getCostMapMarker(visualization_msgs::MarkerArray m, DynamicEDTOctomap *ptr, ros::Publisher pub);
            void get_MMD_Map_Marker( visualization_msgs::MarkerArray m, DynamicEDTOctomap *ptr, ros::Publisher pub );
            void convert_point_to_key(octomap::point3d inpt  ,  int &key_x , int &key_y , int &key_z  );
            std::tuple< int , int , int> convert_point_to_key_external( octomap::point3d inPt);
            double  compute_EDT_interpolation( float distance_at_query_point);
            double compute_EDT_linear_transforms(Eigen::MatrixXf actual_distribution);
            void assign_weights();

            
            Eigen::Vector3d compute_gradients_MMD_map(octomap::point3d  grad_pt );

            struct MMD_Cost{
            float MMD_cost_per_point;
            float mean_obstacle_distance;
            // octomap::point3d  mean_obstacle_location; 
        };

        struct MMD_Gradients{

            float gradient_x;
            float gradient_y;
            float gradient_z; 
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


        struct MMD_grad_key{

            int x_grad_key;
            int y_grad_key;
            int z_grad_key;
            bool operator<(const MMD_grad_key& t2) const
            {
                return (this->x_grad_key < t2.x_grad_key);
            }
        };

        std::map<MMD_Map_key , MMD_Cost> MMD_data;
        std::map<MMD_grad_key , MMD_Gradients> MMD_gradient_data;


        
        private:
            Eigen::Vector3d convertToEigenPt(octomap::point3d pt); // convert octomap::point3d to Eigen::Vector3d
            octomap::point3d convertToOctomapPt(Eigen::Vector3d pt); // convert Eigen::Vector3d to octomap::point3d
            void update_points_for_gradients(octomap::point3d grad_pt);
            void convert_point_to_key_tuple();
            Eigen::MatrixXd  update_coeff();

            std::tuple< int, int, int> point000;
            std::tuple< int, int, int> point100;
            std::tuple< int, int, int> point010;
            std::tuple< int, int, int> point001;

            std::tuple< int, int, int> point110;
            std::tuple< int, int, int> point101;
            std::tuple< int, int, int> point011;

            std::tuple< int, int, int> point111;


            octomap::point3d pt000;
            octomap::point3d pt100;
            octomap::point3d pt010;
            octomap::point3d pt001;

            octomap::point3d pt110;
            octomap::point3d pt011;
            octomap::point3d pt101;

            octomap::point3d pt111;

            Eigen::Vector3d temp_key_point000  ;
            Eigen::Vector3d temp_key_point100  ;
            Eigen::Vector3d temp_key_point010  ;
            Eigen::Vector3d temp_key_point001  ;
            Eigen::Vector3d temp_key_point110  ;
            Eigen::Vector3d temp_key_point011  ;
            Eigen::Vector3d temp_key_point101  ;
            Eigen::Vector3d temp_key_point111  ;                                




 
    };
}




std::tuple< int , int , int> Map3D::OctoMapEDT::convert_point_to_key_external( octomap::point3d inPt) 
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


void  Map3D::OctoMapEDT::convert_point_to_key( octomap::point3d  inpt , int &key_x , int &key_y , int &key_z )
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



//////////////////////////////////////////////////////////////////////////////////////
/** set start and end of Map to (0,0,0) **/
Map3D::OctoMapEDT::OctoMapEDT()
{
    std::cout<<"_____Mapping initialized ..."<<std::endl;

    std::cout<<"\n";

    edtMapRange =  5.0; 
    
    mapUpdateThreshold = 0.25;


}


///////////////////////////////////////////////////////////////////////////////////////
/** delete the trees to save memory **/
void Map3D::OctoMapEDT::erase()
{
    delete tree;
    delete new_tree;
}


/////////////////////////////////////////////////////////////////////////////////////
/** get minimum and maximum in an octomap **/
void Map3D::OctoMapEDT::setMinMax()
{
    double x,y,z;
    tree->getMetricMax(x,y,z);

    max.x() = x;
    max.y() = y;
    max.z() = z; 

    tree->getMetricMin(x,y,z);
    min.x() = x;
    min.y() = y;
    min.z() = z;
}


/////////////////////////////////////////////////////////////////////////////////////
/** Check if the drone has travelled enough distance to update the map again **/
bool Map3D::OctoMapEDT::ifUpdateMap(Eigen::Vector3d pt1)
{
    Eigen::Vector3d pt = convertToEigenPt(mapStart);
    Eigen::Vector3d diff = pt1 - pt;
    
    if(diff.norm()>mapUpdateThreshold)
    {
        std::cout<<">>>>>>>>>>>> Map Updated from location "<<pt1.transpose()<< "to "<<pt.transpose()<<std::endl;
        std::cout<<"\n";
        return true;
    }
    else
    {
        return false;
    }
    
}


/////////////////////////////////////////////////////////////////////////////////////
/** set mapStart to the point where mapping starts/re-initializes **/
void Map3D::OctoMapEDT::setStartPosition(Eigen::Vector3d pt)
{
    mapStart = convertToOctomapPt(pt);
}


////////////////////////////////////////////////////////////////////////////////////
/**
 * set the start and end points of the map using the current pose, 
 * max ranges along each axis and min/max of octomap 
 **/

void Map3D::OctoMapEDT::setMapRange(Eigen::Vector3d pt)
{
    octomap::point3d p = convertToOctomapPt(pt);

    std::cout<<"Octomap point "<<p<<std::endl;
    
   /** set start coordinates **/
    start.x() = p.x() - edtMapRange;
    start.y() = p.y() - edtMapRange;;
    start.z() = p.z() - edtMapRange;

   /** set end coordinates **/
    end.x() = p.x() + edtMapRange;
    end.y() = p.y() + edtMapRange;
    end.z() = p.z() + edtMapRange;

    std::cout<<"Map from "<<start<<" to "<<end<<std::endl;

}

/////////////////////////////////////////////////////////////////////////////////////////
/** check if a point lies in map or not **/
bool Map3D::OctoMapEDT::isInMap(octomap::point3d pt)
{
    if ( start.x() <= pt.x() <= end.x() && start.y() <= pt.y() <= end.y() && start.z() <= pt.z() <= end.z())
        return true;
    else
        return false;
}


//////////////////////////////////////////////////////////////////////////////////////////
/** publish the costmap obtained from the octree **/
void Map3D::OctoMapEDT::getCostMapMarker(visualization_msgs::MarkerArray m, DynamicEDTOctomap* ptr, ros::Publisher pub)
{
    unsigned char maxDepth = 16;
    uint32_t shape = visualization_msgs::Marker::CUBE;
    int count  = 0;


    // set a bounding box using the start and end variables
    for(octomap::OcTree::leaf_bbx_iterator it = tree->begin_leafs_bbx(start, end, maxDepth), bbx_end = tree->end_leafs_bbx(); it != bbx_end; std::advance(it,3))
    {
        // std::cout<<it.getCoordinate()<<std::endl;
        octomap::point3d pt = it.getCoordinate();

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
                        
        marker.color.r = (1 - dist/10.0)*(1 - dist/10.0);
        marker.color.g = std::sqrt(std::sqrt(dist/10.0));
        marker.color.b = 0.0;
        marker.color.a = 0.5;
                    
        marker.lifetime = ros::Duration();
        m.markers.push_back(marker);
        
    }
        pub.publish(m);
}



double get_EDT_cost( Eigen::MatrixXf actual_distribution){

/*

    std::random_device rd{};
    std::mt19937 gen{rd()};

    std::normal_distribution<> d{distance,2};
    int num_samples_of_distance_distribution = 100;

    Eigen::MatrixXf d_dist(1, num_samples_of_distance_distribution);
    Eigen::MatrixXf actual_distribution(1,num_samples_of_distance_distribution);
    Eigen::MatrixXf actual_distribution_transpose(num_samples_of_distance_distribution,1);
    Eigen::MatrixXf actual_distribution_square(1,num_samples_of_distance_distribution);
    Eigen::MatrixXf actual_distribution_square_transpose(num_samples_of_distance_distribution , 1);

    Eigen::MatrixXf desired_distribution(1, num_samples_of_distance_distribution);
    Eigen::MatrixXf alpha_weights(1, num_samples_of_distance_distribution);
    float safe_radius , temp_val1, temp_val2 ;
    safe_radius = 0.75 ;
    double res_val; 

    

    


    for (int i=0 ; i<num_samples_of_distance_distribution; i++){

      temp_val1 = safe_radius -float(d(gen));
      temp_val2 = std::max( float(0.0), temp_val1);
      actual_distribution(0, i) = temp_val2;
      actual_distribution_transpose(i ,0) = temp_val2;

      actual_distribution_square(0,i) = temp_val2*temp_val2;
      actual_distribution_square_transpose(i ,0) = temp_val2*temp_val2;

    //  cout << "till here 0.3" << endl;
      desired_distribution(0, i) = float(0.0); 
      alpha_weights(0, i) = float(1.0);
    }

    Eigen::MatrixXf kernel_matrix1(num_samples_of_distance_distribution, num_samples_of_distance_distribution);
    Eigen::MatrixXf temp_kernel_matrix1(num_samples_of_distance_distribution, num_samples_of_distance_distribution);
    Eigen::MatrixXf temp_kernel_matrix2(num_samples_of_distance_distribution, num_samples_of_distance_distribution);


    Eigen::MatrixXf One_matrix( num_samples_of_distance_distribution, num_samples_of_distance_distribution);
    One_matrix.setOnes();


    temp_kernel_matrix1 = actual_distribution_transpose*actual_distribution;
    temp_kernel_matrix2 = actual_distribution_square_transpose*actual_distribution_square;
    kernel_matrix1 = One_matrix + 2*temp_kernel_matrix1 + temp_kernel_matrix2;

   Eigen::Matrix<float, 1,1>  res; 
    Eigen::Matrix<float, 1,1>  res1; 
    Eigen::Matrix<float, 1,1>  res2; 
    Eigen::Matrix<float, 1,1>  res3; 

    float temp_value =0; 

     


res1 = alpha_weights*kernel_matrix1*(alpha_weights.transpose()) ; 


res2(0 , 0 ) = num_samples_of_distance_distribution;
res3(0,0) = num_samples_of_distance_distribution ;


    res(0,0) = res1(0,0) - 2*res2(0,0) +res3(0,0) ;
    res_val = double( res(0,0) );


    return res_val; 



*/
int num_samples_of_distance_distribution ;
num_samples_of_distance_distribution= actual_distribution.cols();
Eigen::MatrixXf squared_dist(1, num_samples_of_distance_distribution);
Eigen::MatrixXf temp(1, num_samples_of_distance_distribution);
Eigen::MatrixXf alpha_weights(1, num_samples_of_distance_distribution);

alpha_weights.setOnes();


// temp = actual_distribution.array();
squared_dist = (actual_distribution.array()).square();

Eigen::MatrixXf kernel_matrix1(num_samples_of_distance_distribution, num_samples_of_distance_distribution);
Eigen::MatrixXf temp_kernel_matrix1(num_samples_of_distance_distribution, num_samples_of_distance_distribution);
Eigen::MatrixXf temp_kernel_matrix2(num_samples_of_distance_distribution, num_samples_of_distance_distribution);


Eigen::MatrixXf One_matrix( num_samples_of_distance_distribution, num_samples_of_distance_distribution);
One_matrix.setOnes();

// std::cout << " Entered the main funciton /////////////" << std::endl;
temp_kernel_matrix1 = actual_distribution.transpose()*actual_distribution;
temp_kernel_matrix2 = squared_dist.transpose()*squared_dist;
kernel_matrix1 = One_matrix + 2*temp_kernel_matrix1 + temp_kernel_matrix2;

Eigen::Matrix<float, 1,1>  res; 
Eigen::Matrix<float, 1,1>  res1; 
Eigen::Matrix<float, 1,1>  res2; 
Eigen::Matrix<float, 1,1>  res3; 

float temp_value =0; 

     

double res_val;
res1 = alpha_weights*kernel_matrix1*(alpha_weights.transpose()) ; 

res2(0 , 0 ) = num_samples_of_distance_distribution*num_samples_of_distance_distribution;
res3(0,0) = num_samples_of_distance_distribution*num_samples_of_distance_distribution ;

res(0,0) = res1(0,0) - 2*res2(0,0) +res3(0,0) ;
res_val = double( res(0,0) );

return res_val; 


} 




double Map3D::OctoMapEDT::compute_EDT_interpolation( float distance_at_query_point )
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


void Map3D::OctoMapEDT::assign_weights()
{
int rows = 100;
int cols =5;
 std::string file = "../src/CCO_VOXEL/include/CCO_VOXEL/weight.csv" ;
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



double Map3D::OctoMapEDT::compute_EDT_linear_transforms(Eigen::MatrixXf actual_distribution)
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

void Map3D::OctoMapEDT::get_MMD_Map_Marker( visualization_msgs::MarkerArray mdd_marker, DynamicEDTOctomap* ptr, ros::Publisher MMD_map_pub )
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
     assign_weights();
     int mmd_num_counter = 0 ;


    for (int i=0 ; i<num_samples_of_distance_distribution; i++){


      noise_distribution(0 ,i) = radius(0 ,i) - float(noise(gen));
    }



    for(octomap::OcTree::leaf_bbx_iterator it = tree->begin_leafs_bbx(start, end, maxDepth), bbx_end = tree->end_leafs_bbx(); it != bbx_end; std::advance(it,3))
    {




        // std::cout<<it.getCoordinate()<<std::endl;
        octomap::point3d pt = it.getCoordinate();
        // Eigen::Matrix<int , 3,1,> temp_key_vector;

        convert_point_to_key( pt  , key_x , key_y , key_z);

        MMD_Map_key Key_data;

        Key_data.xkey= key_x;
        Key_data.ykey =key_y;
        Key_data.zkey=key_z;

        // MMD_costmap.xkey_temp = key_x;
        // MMD_costmap.ykey_temp = key_y;
        // MMD_costmap.zkey_temp = key_z;




        key_from_octomap = it.getKey();  // conversion from MMD_Map_key to key_from_octomap is key_from_octomap = MMD_Map_key + max_val_tree the mac_val_tree is a constant

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

            MMD_val = compute_EDT_linear_transforms(actual_distribution ); //get_EDT_cost(actual_distribution );  //compute_EDT_interpolation(dist ); //get_EDT_cost(dist ); 

            mmd_num_counter++;
            std::cout <<  mmd_num_counter << std::endl;
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

             // MMD_costmap.MMD_cost_per_point_temp = MMD_val;
             // MMD_costmap.mean_obstacle_distance_temp = dist;

             // MMD_costmap.update_mmd_map();




            // std::cout << MMD_data[Key_data].MMD_cost_per_point <<  " " << MMD_val << std::endl;

            // std::cout << Key_data.xkey << " " << Key_data.ykey << "  " << Key_data.zkey << std::endl;

             // std::cout << key_from_octomap[0] << "  "  <<  key_from_octomap[1] << " "<<  key_from_octomap[2] <<  std::endl;
              std::cout << "--------------" << std::endl;


        marker.color.r =  std::sqrt(std::sqrt(MMD_val/1000.0));
        marker.color.g = (1 - MMD_val/1000.0)*(1 - MMD_val/1000.0);
        marker.color.b = 0.0;
        marker.color.a = (1 - dist/10.0)*(1 - dist/10.0);
                    
        marker.lifetime = ros::Duration();
        mdd_marker.markers.push_back(marker);



    }

    MMD_map_pub.publish(mdd_marker);
    delta_Vector.x() = gradient_stepx;
    delta_Vector.y() = gradient_stepy;
    delta_Vector.z() = gradient_stepz;
    delta = convertToOctomapPt(delta_Vector);
    Eigen::Vector3d gradient_Vector_per_point;
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

Eigen::Vector3d Map3D::OctoMapEDT::compute_gradients_MMD_map(octomap::point3d  grad_pt )
{

    Eigen::MatrixXd point_coeff(8,1);

    update_points_for_gradients(grad_pt); 
    point_coeff = update_coeff();
    double gradx, grady, gradz ;
    Eigen::Vector3d gradient_Vector;

    

    gradx = point_coeff(1,0) + point_coeff(4,0)*pt000.y()  + point_coeff(5,0)*pt000.z() +point_coeff(7,0)*pt000.z()*pt000.y();
    grady = point_coeff(2,0) + point_coeff(4,0)*pt000.x()  + point_coeff(6,0)*pt000.z() +point_coeff(7,0)*pt000.x()*pt000.z();
    gradz = point_coeff(3,0) + point_coeff(5,0)*pt000.x()  + point_coeff(6,0)*pt000.y() +point_coeff(7,0)*pt000.x()*pt000.y();

    gradient_Vector.x() = gradx;
    gradient_Vector.y() = grady;
    gradient_Vector.z() = gradz;

    return gradient_Vector;
}


Eigen::MatrixXd Map3D::OctoMapEDT::update_coeff()
{
    convert_point_to_key_tuple();

    Eigen::MatrixXd coeff(8,1);
    Eigen::MatrixXd A(8,8);
    Eigen::MatrixXd B(8,1);

    MMD_Map_key MMK;

    for(int i =0 ; i < 8 ;i++){

            A(i ,0) = 1;         
    }


    A(0,1) = temp_key_point000.x();
    A(0,2) = temp_key_point000.y();
    A(0,3) = temp_key_point000.z();

    A(0,4) = temp_key_point000.x()*temp_key_point000.y();
    A(0,5 ) = temp_key_point000.x()*temp_key_point000.z() ;

    A(0,6 ) = temp_key_point000.y()*temp_key_point000.z() ;
    A(0,7 ) = temp_key_point000.x()*temp_key_point000.y()*temp_key_point000.z() ;



    A(1,1) = temp_key_point100.x();
    A(1,2) = temp_key_point100.y();
    A(1,3) = temp_key_point100.z();

    A(1,4) = temp_key_point100.x()*temp_key_point100.y();
    A(1,5 ) = temp_key_point100.x()*temp_key_point100.z() ;

    A(1,6 ) = temp_key_point100.y()*temp_key_point100.z() ;
    A(1,7 ) = temp_key_point100.x()*temp_key_point100.y()*temp_key_point100.z() ;




    A(2,1) = temp_key_point010.x();
    A(2,2) = temp_key_point010.y();
    A(2,3) = temp_key_point010.z();

    A(2,4) = temp_key_point010.x()*temp_key_point010.y();
    A(2,5 ) = temp_key_point010.x()*temp_key_point010.z() ;

    A(2,6 ) = temp_key_point010.y()*temp_key_point010.z() ;
    A(2,7 ) = temp_key_point010.x()*temp_key_point010.y()*temp_key_point010.z() ;


    A(3,1) = temp_key_point010.x();
    A(3,2) = temp_key_point010.y();
    A(3,3) = temp_key_point010.z();

    A(3,4) = temp_key_point010.x()*temp_key_point010.y();
    A(3,5 ) = temp_key_point010.x()*temp_key_point010.z() ;

    A(3,6 ) = temp_key_point010.y()*temp_key_point010.z() ;
    A(3,7 ) = temp_key_point010.x()*temp_key_point010.y()*temp_key_point010.z() ;


    A(4,1) = temp_key_point110.x();
    A(4,2) = temp_key_point110.y();
    A(4,3) = temp_key_point110.z();

    A(4,4) = temp_key_point110.x()*temp_key_point110.y();
    A(4,5 ) = temp_key_point110.x()*temp_key_point110.z() ;

    A(4,6 ) = temp_key_point110.y()*temp_key_point110.z() ;
    A(4,7 ) = temp_key_point110.x()*temp_key_point110.y()*temp_key_point110.z() ;



    A(5,1) = temp_key_point101.x();
    A(5,2) = temp_key_point101.y();
    A(5,3) = temp_key_point101.z();

    A(5,4) = temp_key_point101.x()*temp_key_point101.y();
    A(5,5 ) = temp_key_point101.x()*temp_key_point101.z() ;

    A(5,6 ) = temp_key_point101.y()*temp_key_point101.z() ;
    A(5,7 ) = temp_key_point101.x()*temp_key_point101.y()*temp_key_point101.z() ;


    A(6,1) = temp_key_point011.x();
    A(6,2) = temp_key_point011.y();
    A(6,3) = temp_key_point011.z();

    A(6,4) = temp_key_point011.x()*temp_key_point011.y();
    A(6,5 ) = temp_key_point011.x()*temp_key_point011.z() ;

    A(6,6 ) = temp_key_point011.y()*temp_key_point011.z() ;
    A(6,7 ) = temp_key_point011.x()*temp_key_point011.y()*temp_key_point011.z() ;

    A(7,1) = temp_key_point111.x();
    A(7,2) = temp_key_point111.y();
    A(7,3) = temp_key_point111.z();

    A(7,4) = temp_key_point111.x()*temp_key_point111.y();
    A(7,5 ) = temp_key_point111.x()*temp_key_point111.z() ;

    A(7,6 ) = temp_key_point111.y()*temp_key_point111.z() ;
    A(7,7 ) = temp_key_point111.x()*temp_key_point111.y()*temp_key_point111.z() ;


    MMK.xkey = int(std::get<0>(point000));
    MMK.ykey = int(std::get<1>(point000));
    MMK.zkey = int(std::get<2>(point000));
    B(0,0) = double(MMD_data[MMK].MMD_cost_per_point);


    MMK.xkey = int(std::get<0>(point100));
    MMK.ykey = int(std::get<1>(point100));
    MMK.zkey = int(std::get<2>(point100));

    B(1,0) = double(MMD_data[MMK].MMD_cost_per_point);


    MMK.xkey = int(std::get<0>(point010));
    MMK.ykey = int(std::get<1>(point010));
    MMK.zkey = int(std::get<2>(point010));
    B(2,0) = double(MMD_data[MMK].MMD_cost_per_point);


    MMK.xkey = int(std::get<0>(point001));
    MMK.ykey = int(std::get<1>(point001));
    MMK.zkey = int(std::get<2>(point001));
    B(3,0) = double(MMD_data[MMK].MMD_cost_per_point);



    MMK.xkey = int(std::get<0>(point110));
    MMK.ykey = int(std::get<1>(point110));
    MMK.zkey = int(std::get<2>(point110));
    B(4,0) = double(MMD_data[MMK].MMD_cost_per_point);



    MMK.xkey = int(std::get<0>(point101));
    MMK.ykey = int(std::get<1>(point101));
    MMK.zkey = int(std::get<2>(point101));
    B(5,0) = double(MMD_data[MMK].MMD_cost_per_point);


    MMK.xkey = int(std::get<0>(point011));
    MMK.ykey = int(std::get<1>(point011));
    MMK.zkey = int(std::get<2>(point011));
    B(6,0) = double(MMD_data[MMK].MMD_cost_per_point);


    MMK.xkey = int(std::get<0>(point111));
    MMK.ykey = int(std::get<1>(point111));
    MMK.zkey = int(std::get<2>(point111));
    B(7,0) = double(MMD_data[MMK].MMD_cost_per_point);



// coeff = A.colPivHouseholderQr().solve(B);

    coeff = A.completeOrthogonalDecomposition().solve(B);


return coeff;

}



void  Map3D::OctoMapEDT::convert_point_to_key_tuple()
{


    temp_key_point000 = convertToEigenPt( pt000);
    temp_key_point100 = convertToEigenPt(pt100);
    temp_key_point010  = convertToEigenPt(pt010);
    temp_key_point001  = convertToEigenPt(pt001);

     temp_key_point110 = convertToEigenPt( pt110);
    temp_key_point011 = convertToEigenPt(pt011);
    temp_key_point101  = convertToEigenPt(pt101);

    temp_key_point111 = convertToEigenPt( pt111);



    double x000 , y000 ,z000;


    double x100 , y100 ,z100;
    double x010 , y010 ,z010;
    double x001 , y001 ,z001;

    double x110 , y110 ,z110;
    double x101 , y101 ,z101;
    double x011 , y011 ,z011;

    double x111 , y111 ,z111;


    x000 = double( temp_key_point000.x() );
    y000 = double ( temp_key_point000.y());
    z000 = double( temp_key_point000.z());


    x100 = double( temp_key_point100.x() );
    y100 = double ( temp_key_point100.y());
    z100 = double( temp_key_point100.z());


    x010 = double( temp_key_point010.x() );
    y010 = double ( temp_key_point010.y());
    z010 = double( temp_key_point010.z());

    x001 = double( temp_key_point001.x() );
    y001 = double ( temp_key_point001.y());
    z001 = double( temp_key_point001.z());

    x110 = double( temp_key_point110.x() );
    y110 = double ( temp_key_point110.y());
    z110 = double( temp_key_point110.z());

    x101 = double( temp_key_point101.x() );
    y101 = double ( temp_key_point101.y());
    z101 = double( temp_key_point101.z());

    x011 = double( temp_key_point011.x() );
    y011 = double ( temp_key_point011.y());
    z011 = double( temp_key_point011.z());

    x111 = double( temp_key_point111.x() );
    y111 = double ( temp_key_point111.y());
    z111 = double( temp_key_point111.z());


    point000 = std::make_tuple(   (int)floor(x000*resolution_factor) , (int)floor(y000*resolution_factor) , (int)floor(z000*resolution_factor ) ) ;
    point100 = std::make_tuple(   (int)floor(x100*resolution_factor) , (int)floor(y100*resolution_factor) , (int)floor(z100*resolution_factor ) ) ;
    point010 = std::make_tuple(   (int)floor(x010*resolution_factor) , (int)floor(y010*resolution_factor) , (int)floor(z010*resolution_factor ) ) ;
    point001 = std::make_tuple(   (int)floor(x001*resolution_factor) , (int)floor(y001*resolution_factor) , (int)floor(z001*resolution_factor ) ) ;

    point110 = std::make_tuple(   (int)floor(x110*resolution_factor) , (int)floor(y110*resolution_factor) , (int)floor(z110*resolution_factor ) ) ;
    point011 = std::make_tuple(   (int)floor(x101*resolution_factor) , (int)floor(y101*resolution_factor) , (int)floor(z101*resolution_factor ) ) ;
    point101 = std::make_tuple(   (int)floor(x011*resolution_factor) , (int)floor(y011*resolution_factor) , (int)floor(z011*resolution_factor ) ) ;

    point111 = std::make_tuple(   (int)floor(x111*resolution_factor) , (int)floor(y111*resolution_factor) , (int)floor(z111*resolution_factor ) ) ;
}



void Map3D::OctoMapEDT::update_points_for_gradients(octomap::point3d grad_pt)
{

pt000.x() = grad_pt.x();
pt000.y() = grad_pt.y();
pt000.z() = grad_pt.z();

pt100.x() = grad_pt.x() + delta.x();
pt100.y() = grad_pt.y();
pt100.z() = grad_pt.z();


pt010.x() = grad_pt.x();
pt010.y() = grad_pt.y() + delta.y();
pt010.z() = grad_pt.z();

pt001.x() = grad_pt.x();
pt001.y() = grad_pt.y();
pt001.z() = grad_pt.z() +delta.z();


pt110.x() = grad_pt.x() + delta.x();
pt110.y() = grad_pt.y() + delta.y();
pt110.z() = grad_pt.z();


pt101.x() = grad_pt.x() + delta.x();
pt101.y() = grad_pt.y();
pt101.z() = grad_pt.z() + delta.z();


pt011.x() = grad_pt.x();
pt011.y() = grad_pt.y() + delta.y();
pt011.z() = grad_pt.z() + delta.z();


pt100.x() = grad_pt.x() + delta.x();
pt100.y() = grad_pt.y() + delta.y();
pt100.z() = grad_pt.z()+ delta.z();

}

//////////////////////////////////////////////////////////////////////////////////////////
/** convert to octomap::point3d from Eigen::Vector3d **/
octomap::point3d Map3D::OctoMapEDT::convertToOctomapPt(Eigen::Vector3d pt)
{
    octomap::point3d p(pt(0), pt(1), pt(2));
    return p;
}


////////////////////////////////////////////////////////////////////////////////////////
/** convert to Eigen::Vector3d from octomap::point3d **/
Eigen::Vector3d Map3D::OctoMapEDT::convertToEigenPt(octomap::point3d pt)
{
    Eigen::Vector3d p(pt.x(), pt.y(), pt.z());
    return p;
}
