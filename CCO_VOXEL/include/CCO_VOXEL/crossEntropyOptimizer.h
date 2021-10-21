/******** Cross entropy optimization ******/

/************************************************************************
* Take trajectory coeffs and matrices
* Perturb the coefficients
* Check the MMD cost of each trajectory and select the top best of those 
*************************************************************************/

#include"utils.h"
#include"bernstein.h"
#include"bsplineNonUnif.h"
#include"Map.h"
#include<random>
#include<algorithm>
#include"visulization.h"
#include <fstream>
#include <filesystem>



Visualizer::Trajectory_visualizer traj_vis;


namespace Optimizer
{
    class CrossEntropyOptimizer
    {
        public:
            int numIterations;
            int topSamples     = 20;
            double safeRadius  = 1.0;
            int ptsPerTraj     = 50;
            int numSampleTrajs = 50;
            double infCost = 1000000.0;
            Eigen::MatrixXd optimTrajCoeffs = Eigen::MatrixXd::Zero(11,3);
            Eigen::Vector3d var_vector;
            std::string path_to_weights2 ; 
            CrossEntropyOptimizer();
            CrossEntropyOptimizer(int numIterations_);
            std::vector<Eigen::Vector3d> optimizeTrajectory(Bernstein::BernsteinPath bTraj, std::vector<Eigen::Vector3d> wayPts, float execTime, Map3D::OctoMapEDT costMap3D , ros::Publisher sample_trajectory_pub , 
                ros::Publisher plan_dur_pub , std::string path_to_weights);
            double costPerTrajectory(std::vector<Eigen::Vector3d> trajectory, std::vector<Eigen::Vector3d> trajectoryAcc, std::vector<Eigen::Vector3d> initTrajectory, Map3D::OctoMapEDT costMap3D , bool is_mean , 
                ros::Publisher plan_dur_pub);
            double get_variance(Eigen::MatrixXd one_dimension_trajectory , int iter);
            double get_acc_cost( Eigen::Vector3d acc_in );
            double get_elastic_band_cost( std::vector<Eigen::Vector3d> traj_in ); 
            double mmdPerPoint_interpolation(float distance);
            float mmdPerPoint_transforms(  Eigen::MatrixXf actual_distribution);
            void assign_weights();
            float MMD_transformed_features_RBF( Eigen::MatrixXf actual_distribution);
            float RBF_kernel( float val1 , float val2);
            double get_total_smoothness_cost( std::vector<Eigen::Vector3d> optimTraj , float execTime);

            Eigen::Matrix< float , 100, 5> Weights;
            std::ofstream dist_measurments ;
        
        private:
            inline double mmdPerPoint(std::vector<double> actualDistribution, std::vector<double> idealDistribution, std::vector<double> weights, int numEdtSamples);
            inline Eigen::MatrixXd convertVecTrajToMatTraj(std::vector<Eigen::Vector3d> arr);
            inline std::vector<Eigen::Vector3d> convertMatTrajToVecTraj(Eigen::MatrixXd mat);
            
    };

}

/*******************************
* Default constructor          *
********************************/
Optimizer::CrossEntropyOptimizer::CrossEntropyOptimizer()
{
    // std::cout<<"Setting default iterations as 3"<<std::endl;
    
     

    numIterations = 5;
}

/*********************************
* Parameterized constructor      *
*********************************/
Optimizer::CrossEntropyOptimizer::CrossEntropyOptimizer(int numIterations_)
{
    std::cout<<"Optimizer Iterations set to "<<numIterations_<<std::endl;
    numIterations = numIterations_;
}


double Optimizer::CrossEntropyOptimizer::get_total_smoothness_cost( std::vector<Eigen::Vector3d> optimTraj , float execTime)
{

float time =0;
int num_waypts = optimTraj.size() ; 
float delta_t = execTime/num_waypts;
double integral_jerk_cost =0 ;
float prev_t = 0 ; 
for (int i =1 ; i< num_waypts ; i++){
time = (float(i)/float(num_waypts-1))*execTime ;  
// std::cout << time <<std::endl;
delta_t = time - prev_t ;
Eigen::Vector3d pt_prev;
Eigen::Vector3d pt_now;

pt_prev =  optimTraj.at(i-1);
pt_now = optimTraj.at(i);

float length =  (pt_prev -pt_now ).norm() ; 

integral_jerk_cost += std::pow( length , 2) / ( std::pow (time , 5 ));
prev_t = time ;

}


integral_jerk_cost = integral_jerk_cost;



return integral_jerk_cost ;


}


/*****************************
* main optimizer function    *
******************************/
std::vector<Eigen::Vector3d> Optimizer::CrossEntropyOptimizer::optimizeTrajectory(Bernstein::BernsteinPath bTraj, std::vector<Eigen::Vector3d> wayPts, float execTime, Map3D::OctoMapEDT costMap3D , ros::Publisher sample_trajectory_pub ,
    ros::Publisher plan_dur_pub , std::string path_to_weights)
{
    // generate the initial set of coefficients
    std::cout<<"----Generating bernstein trajectory for "<<wayPts.size()<<" points"<<std::endl;
    std::vector<Eigen::Vector3d> prev_mean_bernstein_trajectory; 
    bTraj.generateCoeffMatrices(wayPts.size(), execTime);
    bTraj.generateTrajCoeffs(wayPts);                                            // this would initialize the trajectory coefficients with those of the fast planner
    bTraj.generateCoeffMatrices(ptsPerTraj, execTime);

    Eigen::MatrixXd initCoeff = convertVecTrajToMatTraj(bTraj.coeffs);           // this takes in a vector of 11 indices and returns a matrix of size ptsPerTrajx3

    Eigen::MatrixXd initWayPts = (bTraj.P)*initCoeff;                            // this returns a 50x3 matrix

    std::vector<Eigen::Vector3d> initBernsteinTraj = convertMatTrajToVecTraj(initWayPts); // returns initial trajectory 

    std::vector<Eigen::Vector3d> coeffs_ = bTraj.coeffs;                         // initial coefficients

    path_to_weights2 = path_to_weights ; 

    assign_weights();

    double var  = 5;
    var_vector.x() = 7;
    var_vector.y() = 7;
    var_vector.z() = 7;

    int num_prev_top_traj = 0.2*topSamples; 

    Eigen::MatrixXd prev_TopX( num_prev_top_traj, ptsPerTraj  );
    Eigen::MatrixXd prev_TopY( num_prev_top_traj , ptsPerTraj  );
    Eigen::MatrixXd prev_TopZ( num_prev_top_traj , ptsPerTraj  );


    

    // steps -> randomly perturb -> generate path -> check for mmd cost -> select the best -> update mean and variance -> recompute the best one
    for(int iter = 0; iter<numIterations; iter++)
    {
        std::vector<std::vector<Eigen::Vector3d> > trajs;
        std::vector<std::vector<Eigen::Vector3d> > trajsAcc;

        std::cout<<"Cross entropy Iteration "<<iter<<std::endl;
        
        // perturb the coefficients now

        // std::cout << coeffs_.size() <<  "*************************** "  << std::endl; 
        std::vector<Eigen::MatrixXd> perturbedCoeffs = bTraj.generatePerturbedCoeffs(numSampleTrajs, coeffs_, var_vector);

        
        // std::cout<<"\n *************************************************************************************** \n"<< " done perturbing"  << std::endl;

        // generate the trajectories (xPts, yPts, zPts each of size -> numSampleTrajs x ptsPerTraj)
        Eigen::MatrixXd xPts = ((bTraj.P)*(perturbedCoeffs.at(0).transpose())).transpose();
        Eigen::MatrixXd yPts = ((bTraj.P)*(perturbedCoeffs.at(1).transpose())).transpose();
        Eigen::MatrixXd zPts = ((bTraj.P)*(perturbedCoeffs.at(2).transpose())).transpose();


        // std::cout<<"\n *************************************************************************************** \n"<< " done perturbing 2"  << std::endl;



        /**
        *   also compute stability costs -> by computing jerk and snaps
        *   minimize the snap over time throughout the trajectory
        *   smoothen the trajectory over time
        *   re-allocate time if the jerk values are high
        **/

        Eigen::MatrixXd xAccPts = ((bTraj.Pddot)*(perturbedCoeffs.at(0).transpose())).transpose();
        Eigen::MatrixXd yAccPts = ((bTraj.Pddot)*(perturbedCoeffs.at(1).transpose())).transpose();
        Eigen::MatrixXd zAccPts = ((bTraj.Pddot)*(perturbedCoeffs.at(2).transpose())).transpose();

         // std::cout<<"\n *************************************************************************************** \n"<< " done perturbing 3"  << std::endl;
       
        std::vector<double> costTrajs(numSampleTrajs);

        for(int i = 0; i<numSampleTrajs; i++)
        {
            std::vector<Eigen::Vector3d> traj;
            std::vector<Eigen::Vector3d> trajAcc;
            bool getCost = true;

            for(int j = 0; j<ptsPerTraj; j++)
            {   
                Eigen::Vector3d pt(xPts(i,j), yPts(i,j), zPts(i,j));
                Eigen::Vector3d ptAcc(xAccPts(i,j), yAccPts(i,j), zAccPts(i,j));

                octomap::point3d octoPt(pt(0), pt(1), pt(2));

                if(costMap3D.isInMap(octoPt))
                {
                    float dist_ = costMap3D.costMap->getDistance(octoPt);

                    if(dist_ < 0)
                    {
                        costTrajs.at(i) += infCost;
                        getCost = false;
                    }
                }
                
                if(!costMap3D.isInMap(octoPt))
                {
                    costTrajs.at(i) += infCost;
                    getCost = false;
                }

                traj.push_back(pt);
                trajAcc.push_back(ptAcc);
            }

            trajs.push_back(traj);
            trajsAcc.push_back(trajAcc);

            // now compute cost for each trajectory which is in the map and whose even 1 point does not collide with obstacles
            if(getCost)
            {
                bool is_mean = false;
                double cost = costPerTrajectory(traj, trajAcc, initBernsteinTraj, costMap3D , is_mean , plan_dur_pub);
                costTrajs.at(i) = cost;
            }
        
        }

    std::vector<double> costTrajsorted = costTrajs;
    std::sort(costTrajsorted.begin(), costTrajsorted.end());

    double valOptim = costTrajsorted.at(0);
    auto itrOptim   = std::find(costTrajs.begin(), costTrajs.end(), valOptim);
    int indexOptim  = itrOptim - costTrajs.begin();





    for(int k = 0; k<11; k++)
    {
        optimTrajCoeffs(k,0) = perturbedCoeffs.at(0)(indexOptim,k);//coeffs_.at(k)(0);//perturbedCoeffs.at(0)(indexOptim,k);// // // 
        optimTrajCoeffs(k,1) = perturbedCoeffs.at(1)(indexOptim,k);//coeffs_.at(k)(1);//perturbedCoeffs.at(1)(indexOptim,k);// // //;
        optimTrajCoeffs(k,2) = perturbedCoeffs.at(2)(indexOptim,k);//coeffs_.at(k)(2);//perturbedCoeffs.at(2)(indexOptim,k);// ///
    }

    // std::cout<<"\n"<<std::endl;

    std::vector<Eigen::Vector3d> newCoeffs(11);
    
    std::vector<int> topIndexes;


    if( iter > 0)
    {
    

    for(int index_top = 0; index_top < num_prev_top_traj ; index_top++ ){
 
       xPts.row(index_top)  = prev_TopX.row(index_top);
       yPts.row(index_top)  = prev_TopY.row(index_top);
       zPts.row(index_top)  = prev_TopZ.row(index_top);


    }
}

    for(int i = 0; i<topSamples; i++)
    {
        auto it = costTrajsorted.begin();
        double costVal = *it;
        auto itr = std::find(costTrajs.begin(), costTrajs.end(), costVal);
        int index = itr - costTrajs.begin();
        topIndexes.push_back(index);

        costTrajsorted.erase(costTrajsorted.begin());
    }

    // add up the matrices which are in the topSamples
    Eigen::MatrixXd sum_top_X = Eigen::MatrixXd::Zero(1,11);
    Eigen::MatrixXd sum_top_Y = Eigen::MatrixXd::Zero(1,11);
    Eigen::MatrixXd sum_top_Z = Eigen::MatrixXd::Zero(1,11);

    Eigen::MatrixXd diffTopCoeffsX = Eigen::MatrixXd::Zero(1,11);
    Eigen::MatrixXd diffTopCoeffsY = Eigen::MatrixXd::Zero(1,11);
    Eigen::MatrixXd diffTopCoeffsZ = Eigen::MatrixXd::Zero(1,11);


    Eigen::MatrixXd xCoeff = perturbedCoeffs.at(0);
    Eigen::MatrixXd yCoeff = perturbedCoeffs.at(1);
    Eigen::MatrixXd zCoeff = perturbedCoeffs.at(2);



    Eigen::MatrixXd TopX( topSamples , ptsPerTraj  );
    Eigen::MatrixXd TopY( topSamples , ptsPerTraj  );
    Eigen::MatrixXd TopZ( topSamples , ptsPerTraj  );

    Eigen::MatrixXd TopX_mean( 1 , ptsPerTraj  );
    Eigen::MatrixXd TopY_mean( 1 , ptsPerTraj  );
    Eigen::MatrixXd TopZ_mean( 1 , ptsPerTraj  );


    int index = 0;


     for( int p =0 ; p < topSamples ; p ++)
    {
        TopX.row(p)= xPts.row( topIndexes.at(p));  // xPts dimension numsamples x pointspertraj 
        TopY.row(p)= yPts.row( topIndexes.at(p));
        TopZ.row(p)= zPts.row( topIndexes.at(p));

        if( p < num_prev_top_traj){


        prev_TopX.row(p) = TopX.row(p);
        prev_TopY.row(p) = TopY.row(p);
        prev_TopZ.row(p) = TopZ.row(p);

    }
    } // TopX and TopY and TOpZ dimension is topsamples x pointspertraj 




    traj_vis.visulize_sampled_trajectories( TopX , TopY , TopZ , topSamples , ptsPerTraj,  sample_trajectory_pub);

    std::vector<Eigen::Vector3d> mean_bernstein_trajectory; 

    for( int col =0 ; col < ptsPerTraj ; col++)
    {
        Eigen::Vector3d waypoint_mean; 
    for( int row= 0 ; row < topSamples ; row++){

            TopX_mean(0,col) += TopX( row, col); 
            TopY_mean(0,col) += TopY( row, col); 
            TopZ_mean(0,col) += TopZ( row, col); 
        }

        TopX_mean /= topSamples;
        TopY_mean /= topSamples;
        TopZ_mean /= topSamples;


        waypoint_mean(0) =    TopX_mean(0,col);
        waypoint_mean(1) =  TopY_mean(0,col);
        waypoint_mean(2) =   TopZ_mean(0,col);

    
        mean_bernstein_trajectory.push_back(waypoint_mean);

    }



/* ----------------------------------  */


    Eigen::MatrixXd sumTopCoeffsX = Eigen::MatrixXd::Zero(1,11);
    Eigen::MatrixXd sumTopCoeffsY = Eigen::MatrixXd::Zero(1,11);
    Eigen::MatrixXd sumTopCoeffsZ = Eigen::MatrixXd::Zero(1,11);

    Eigen::MatrixXd Coeffx_top_Set =  Eigen::MatrixXd::Zero( topSamples , 11  );
    Eigen::MatrixXd Coeffy_top_Set = Eigen::MatrixXd::Zero( topSamples , 11  );
    Eigen::MatrixXd Coeffz_top_Set = Eigen::MatrixXd::Zero( topSamples , 11  );

    
    for(int i = 0; i<topSamples; i++)
    {   
        // std::cout<<xCoeff.row(0)<<std::endl;
        
        int index_ = topIndexes.at(i);

        sumTopCoeffsX += xCoeff.row(index_); 
        sumTopCoeffsY += yCoeff.row(index_);
        sumTopCoeffsZ += zCoeff.row(index_);


    }

    /*
    std::cout<<"#####################################################################################################"<<std::endl;
    std::cout<<"At iteration "<<iter<<"Sums are "<<sumTopCoeffsX<<"\n"<<sumTopCoeffsY<<"\n"<<sumTopCoeffsZ<<std::endl;
    std::cout<<"#####################################################################################################"<<std::endl;
    */
        // update the coefficient matrix to be used next time
        for(int j = 0; j<11; j++)
        {

            newCoeffs.at(j)(0) = double(sumTopCoeffsX(0,j)/float(topSamples)); 
            newCoeffs.at(j)(1) = double(sumTopCoeffsY(0,j)/float(topSamples)); 
            newCoeffs.at(j)(2) = double(sumTopCoeffsZ(0,j)/float(topSamples));     
        }






/* -------------------------*/ 


        // now update the variance
        
        var_vector.x() = get_variance(TopX , iter ); // this might be the bug the dimension is 20x50 
        var_vector.y() = get_variance(TopY, iter );
        var_vector.z() = get_variance(TopZ , iter);






        std::cout << var_vector.x()  <<   " " << var_vector.y()  << "  " << var_vector.z()   <<  "  " << "updated variance" << std::endl;

        // variance <<  var_vector.x()  <<   " " << var_vector.y()  << "  " << var_vector.z()  << std::endl;

        bTraj.generateTrajCoeffs(mean_bernstein_trajectory);
        coeffs_.clear(); 
        coeffs_ = bTraj.coeffs;
        // bTraj.coeffs = newCoeffs ;
        // coeffs_ = bTraj.coeffs;

        std::vector<Eigen::Vector3d> mean_trajAcc;


        Eigen::Vector3d Pt_Acc; 

        for( int i =0 ; i < ptsPerTraj ; i++){

           Eigen::Vector3d  Pt_Acc( 0 , 0 , 0 );
           mean_trajAcc.push_back(Pt_Acc);
        }

        // std::cout << coeffs_.size()  << "  " << coeffs_.at(0).rows() << "   " << coeffs_.at(0).cols()  << "  " << bestTraj_accx.rows() << " " << bestTraj_accx.cols()  <<  std::endl;
        

        std_msgs::Float64 dummy_data;
        dummy_data.data = 1000.000;

        for( int j =0 ; j <30 ; j++){

        plan_dur_pub.publish(dummy_data);}

        // ros::Duration(3).sleep();
        bool is_mean = true ;
        double mean_traj_cost = costPerTrajectory(mean_bernstein_trajectory , mean_trajAcc , initBernsteinTraj , costMap3D  ,  is_mean , plan_dur_pub);
        // is_mean =false ; 
        // std::cout << " Iteration Complete change data file name " << std::endl; 

         // ros::Duration(30).sleep();
        double mean_smoothness = get_total_smoothness_cost(mean_bernstein_trajectory , execTime );




    }

    Eigen::MatrixXd bestTraj = ((bTraj.P) * optimTrajCoeffs);


    
    // std::cout<<"\n ***********BEST Traj********** \n"<<bestTraj<<"\n ********************** \n"<<std::endl;

    // convert this matrix to a vector
    // std::cout<<"Best trajectory has dimensions "<<bestTraj.rows()<<" "<<bestTraj.cols()<<std::endl;

    std::vector<Eigen::Vector3d> optimTraj;

    for(int ind = 0; ind < bestTraj.rows(); ind++)
    {
        Eigen::Vector3d pt;

        pt(0) = bestTraj(ind, 0);
        pt(1) = bestTraj(ind, 1);
        pt(2) = bestTraj(ind, 2);

        optimTraj.push_back(pt);

        // std::cout<<pt(0)<<"\t"<<pt(1)<<"\t"<<pt(2)<<std::endl;
    }

    double smoothness_cost = get_total_smoothness_cost(optimTraj , execTime);

    // std::cout << smoothness_cost << "smoothness_cost" <<std::endl;

    return optimTraj;
}



double Optimizer::CrossEntropyOptimizer::get_variance( Eigen::MatrixXd one_dimension_trajectory , int iter )
{

int coeff_size = 11;  
double variance_value=0;
Eigen::MatrixXd mean(1,11);
Eigen::MatrixXd mean_sum(1,1);
mean_sum.setZero();
// std::cout << "here " << one_dimension_trajectory.rows() <<  " "  << one_dimension_trajectory.cols()  << std::endl;
 int total_size = one_dimension_trajectory.rows()*one_dimension_trajectory.cols();

 

for(int k=0 ; k< one_dimension_trajectory.cols() ; k++ ){
for(int l=0 ; l< one_dimension_trajectory.rows() ; l++){
mean_sum(0,0) += one_dimension_trajectory( l, k )/(total_size );
}

}

for( int i=0 ;i < one_dimension_trajectory.rows() ; i++){
for( int j =0; j <11 ; j++){
variance_value += pow( ( one_dimension_trajectory( i,j) - mean_sum(0,0) ), 2);
// std::cout << variance_value << std::endl; 

}
}

variance_value = variance_value/(total_size); 

// std::cout << variance_value << std::endl;
return variance_value ;

}


double Optimizer::CrossEntropyOptimizer::get_acc_cost( Eigen::Vector3d acc_in ){

double acc_val = acc_in.norm();
double amin = -1.0 ; 
double amax = 1.5 ;

double acc_cost =0 ;

if( acc_val > amax || acc_val < amin){

acc_cost = ( acc_val - amin)*(acc_val - amax);
}
else{
    acc_cost =0 ;
}

return acc_cost;
}


double Optimizer::CrossEntropyOptimizer::get_elastic_band_cost( std::vector<Eigen::Vector3d> traj_in )
{
 double elastic_cost=0 ; 
for(int i = 1; i<traj_in.size()-1; i++)
    {
        Eigen::Vector3d pt_before = traj_in.at(i-1);
        Eigen::Vector3d pt_current = traj_in.at(i);
        Eigen::Vector3d pt_next = traj_in.at(i+1);

        Eigen::Vector3d diff; 

        diff = pt_before - 2*pt_current + pt_next ; 
        elastic_cost += diff.norm();

}

return elastic_cost;

}





/************************************************************************************
* Get overall costs for each trajectory
* Overall cost includes collision cost using MMD,stability cost and smoothness cost
************************************************************************************/
double Optimizer::CrossEntropyOptimizer::costPerTrajectory(std::vector<Eigen::Vector3d> traj, std::vector<Eigen::Vector3d> trajAcc, std::vector<Eigen::Vector3d> initBernsteinTraj, Map3D::OctoMapEDT costMap3D , bool is_mean ,
    ros::Publisher plan_dur_pub)
{

    double cost = 0.0;
    double collisionCost  = 0.0;
    double stabilityCost  = 0.0;
    double smoothnessCost = 0.0;
    double elastic_band_cost= 0; 
    int number_of_points_in_distribution =100;

    elastic_band_cost = get_elastic_band_cost(traj);

    // std::ofstream dist_measurments("dist_measurments.csv");

    // Eigen::MatrixXd total_dist_per_iteration(traj.size()*number_of_points_in_distribution  ,1 );

    // std::cout << is_mean <<std::endl;



    for(int i = 0; i<traj.size(); i++)
    {   

        Eigen::Vector3d pt = traj.at(i);
        Eigen::Vector3d ptAcc = trajAcc.at(i);
        Eigen::Vector3d ptInit = initBernsteinTraj.at(i);

        stabilityCost  +=   get_acc_cost(ptAcc);  // ptAcc.norm();
        smoothnessCost += (pt - ptInit).norm();

        /** collision cost calculation **/
        octomap::point3d p(pt(0), pt(1), pt(2));
        float dist =  costMap3D.costMap->getDistance(p);
        if(dist < 2.0)
        {
            // generate random distribution around this value
            
            std::default_random_engine de(time(0));
            std::normal_distribution<double> edtDist(dist, 1.0);

            std::vector<double> actualDistibution, idealDistribution, weights;
            Eigen::MatrixXf actual_distribution( 1 , number_of_points_in_distribution);

            for(int r = 0; r<number_of_points_in_distribution; r++)
            {
                actual_distribution(0, r) = (std::max(0.0, (safeRadius - edtDist(de))));
                if(is_mean == true){
                    std_msgs::Float64 distance;
                    distance.data = float(actual_distribution(0, r)) ; 
                    // ros::Duration(0.1).sleep();

                plan_dur_pub.publish(distance) ;
                // std::cout << "mean distance publish" << std::endl;
            }
                
            }

            // collisionCost += mmdPerPoint(actualDistibution, idealDistribution, weights, 50);

            

            // collisionCost += mmdPerPoint_interpolation( dist);



            collisionCost += mmdPerPoint_transforms(actual_distribution); // MMD_transformed_features_RBF



        }
        else
        {   
            collisionCost += 0;

            if(is_mean == true){
                std_msgs::Float64 temp_distance;

                temp_distance.data = 0 ;

                plan_dur_pub.publish(temp_distance) ;

                // std::cout <<  "outside" << std::endl;
            }
        }
    } 

    // dist_measurments.close();

    // std::cout<<"Cost values are: {collision, stability, smoothness} "<<collisionCost<<"\t"<<stabilityCost<<"\t"<<smoothnessCost<<std::endl;
    cost = collisionCost  + 0.50*stabilityCost  + 0.001*elastic_band_cost;
    return cost;
}






float Optimizer::CrossEntropyOptimizer::mmdPerPoint_transforms(  Eigen::MatrixXf actual_distribution)
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

return float(res_val); 


}


void Optimizer::CrossEntropyOptimizer::assign_weights()
{

int rows = 100;
int cols =5;
// std::cout << path_to_weights2 <<" " << " ------------------ " << std::endl;


 std::string file = path_to_weights2 ;
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

  std::cout << Weights.rows() << " " << Weights.cols() << "  " << "Weights assigned" << std::endl;

}





double Optimizer::CrossEntropyOptimizer::mmdPerPoint_interpolation( float distance)
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
        distance_matrix(0 ,i) = pow( distance ,i);
    }

    result = coefficient_matrix*distance_matrix.transpose();

    return double(result(0,0));
}



float Optimizer::CrossEntropyOptimizer::RBF_kernel( float val1, float val2){

float data1 = std::pow( (val1 - val2),  2 )/ (  -2*std::pow( 0.1,  2 ) );
float data2 =  std::exp( data1);

return data2;

}

float Optimizer::CrossEntropyOptimizer::MMD_transformed_features_RBF( Eigen::MatrixXf actual_distribution )
{

Eigen::MatrixXf transformed_features(1,  5);


transformed_features = actual_distribution*Weights;

// std::cout << transformed_features << std::endl;

int num_samples_of_distance_distribution =transformed_features.cols() ;

Eigen::MatrixXf Matrix1( num_samples_of_distance_distribution , num_samples_of_distance_distribution);
Eigen::MatrixXf Matrix2( num_samples_of_distance_distribution , num_samples_of_distance_distribution);
Eigen::MatrixXf Matrix3( num_samples_of_distance_distribution , num_samples_of_distance_distribution);

for( int i =0; i < num_samples_of_distance_distribution ; i++){

  for( int j =0; j < num_samples_of_distance_distribution ;j++){

    Matrix1( i,j )= RBF_kernel( transformed_features( 0, i) , transformed_features(0,j));
    Matrix2( i,j )= RBF_kernel( transformed_features( 0, i) , 0.0);
    Matrix3( i,j )= RBF_kernel( 0.0 ,0.0);
  }


}

Eigen::MatrixXf alpha_weights( 1, num_samples_of_distance_distribution);
alpha_weights.setOnes();



Eigen::Matrix<float, 1,1>  cost1; 
Eigen::Matrix<float, 1,1>  cost2; 
Eigen::Matrix<float, 1,1>  cost3; 


cost1 = alpha_weights*Matrix1*alpha_weights.transpose();
cost2 = alpha_weights*Matrix2*alpha_weights.transpose();
cost3 = alpha_weights*Matrix3*alpha_weights.transpose();

float MMD_RBF_cost = cost1(0,0) - 2*cost2(0,0) + cost3(0,0);

return MMD_RBF_cost;

} 


/*****************************************************
* Get collision cost for each point in the trajectory
******************************************************/
inline double Optimizer::CrossEntropyOptimizer::mmdPerPoint(std::vector<double> actualDistibution, std::vector<double> idealDistribution, std::vector<double> weights, int numEdtSamples)
{
    Eigen::MatrixXd kernel1(numEdtSamples, numEdtSamples);
    Eigen::MatrixXd kernel2(numEdtSamples, numEdtSamples);
    Eigen::MatrixXd kernel3(numEdtSamples, numEdtSamples);
    Eigen::MatrixXd weights_(1,numEdtSamples);

    for(int i = 0; i<numEdtSamples; i++)
    {
        weights_(0,i) = weights.at(i);

        for(int j = 0; j<numEdtSamples; j++)
        {
            kernel1(i,j) = pow((actualDistibution.at(i)*actualDistibution.at(j) + 1), 2);
            kernel2(i,j) = pow((actualDistibution.at(i)*idealDistribution.at(j) + 1), 2);
            kernel3(i,j) = pow((idealDistribution.at(i)*idealDistribution.at(j) + 1), 2);
        }
    }
    Eigen::MatrixXd Res1 = weights_ * kernel1 * (weights_.transpose());
    Eigen::MatrixXd Res2 = weights_ * kernel2 * (weights_.transpose());
    Eigen::MatrixXd Res3 = weights_ * kernel3 * (weights_.transpose());

    return (Res1(0) -2*Res2(0) + Res3(0));
}

/*****************************************************
* Convert vector to array (nx3) size
******************************************************/
inline Eigen::MatrixXd Optimizer::CrossEntropyOptimizer::convertVecTrajToMatTraj(std::vector<Eigen::Vector3d> arr)
{
    Eigen::MatrixXd trajPoints =  Eigen::MatrixXd::Zero(arr.size(), 3);
    
    for(int i = 0; i<arr.size(); i++)
    {
        trajPoints(i,0) = arr.at(i)(0);
        trajPoints(i,1) = arr.at(i)(1);
        trajPoints(i,2) = arr.at(i)(2);
    }

    return trajPoints;
}

/*****************************************************
* Convert (nx3) array to vector
******************************************************/
inline std::vector<Eigen::Vector3d> Optimizer::CrossEntropyOptimizer::convertMatTrajToVecTraj(Eigen::MatrixXd mat)
{
    std::vector<Eigen::Vector3d> trajPoints;

    for(int i = 0; i<mat.rows(); i++)
    {
        Eigen::Vector3d pt;
        pt(0) = mat(i,0);
        pt(1) = mat(i,1);
        pt(2) = mat(i,2);

        trajPoints.push_back(pt);
    }

    return trajPoints;
}

