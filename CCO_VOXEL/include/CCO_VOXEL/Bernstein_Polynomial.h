#include <random>
#include"utils.h"

namespace Bernstein {



	class Bernstein_polynomial_function{

	public:

		int* num_points_in_path;
		int order =10;
		Eigen::MatrixXd *P ; 
		Eigen::MatrixXd *Pdot;
		Eigen::MatrixXd *Pddot;  
		std::vector<Eigen::Vector3d> coeffs;



		void get_10th_order_bernstein_matrix( float execTime , int numWayPts);
		void get_coefficients(Eigen::MatrixXd Path_vector  , int numWayPts );
		// Eigen::MatrixXf get_P_matrix_10th_order(  Eigen::MatrixXf  P,   Eigen::MatrixXf delta_time_vector , int num_waypts);
		// Eigen::MatrixXf get_Pdot_matrix_10th_order(  Eigen::MatrixXf  Pdot,   Eigen::MatrixXf delta_time_vector , int num_waypts);
		// Eigen::MatrixXf get_Pddot_matrix_10th_order(  Eigen::MatrixXf  Pddot,   Eigen::MatrixXf delta_time_vector , int num_waypts);
		float binomialCoefficient(int n, int i);
		float factorial(int x);


	};
}



void Bernstein::Bernstein_polynomial_function::get_coefficients( Eigen::MatrixXd Path_vector , int numWayPts  )
{

// auto i = Bernstein_and_derevative_matrix.begin();
// Eigen::MatrixXf P_in = *i;
// i++;
// Eigen::MatrixXf Pdot_in = *i;
// i++;
// Eigen::MatrixXf Pddot_in = *i;


Eigen::MatrixXd P_in(numWayPts , 11); 
Eigen::MatrixXf Pdot_in(numWayPts , 11); 
Eigen::MatrixXf Pddot_in(numWayPts , 11); 

 P_in= *P;
 Pdot_in = *Pdot;
 Pddot_in = *Pddot;


Eigen::Vector3d X_init; 
X_init = Path_vector.row(0);

Eigen::Vector3d X_fin; 
X_fin = Path_vector.row(numWayPts-1);


Eigen::MatrixXd A_eq(6,11);
A_eq.setZero();
Eigen::MatrixXd b_eq_x(6,1);
Eigen::MatrixXd b_eq_y(6,1);
Eigen::MatrixXd b_eq_z(6,1);

double vx_init = 0.0 ;
double vy_init = 0.0;
double vz_init = 0 ;

double ax_init = 0.0;
double ay_init = 0.0;
double az_init = 0 ;

double vx_fin = 0.0;
double vy_fin = 0.0;
double vz_fin = 0.0;

double ax_fin = 0.0;
double ay_fin = 0.0;
double az_fin = 0.0;

A_eq.row(0) = P_in.row(0);
A_eq.row(1) = Pdot_in.row(0);
A_eq.row(2) = Pddot_in.row(0);
A_eq.row(3) = P_in.row(numWayPts-1);
A_eq.row(4) = Pdot_in.row(numWayPts-1);
A_eq.row(5) = Pddot_in.row(numWayPts-1);



b_eq_x(0,0) = X_init(0);
b_eq_x(1,0) = vx_init;
b_eq_x(2,0) = ax_init;
b_eq_x(3,0) = X_fin(0);
b_eq_x(4,0) = vy_fin;
b_eq_x(5,0) = ay_fin;

b_eq_y(0,0) = X_init(1);
b_eq_y(1,0) = vy_init;
b_eq_y(2,0) = ay_init;
b_eq_y(3,0) = X_fin(1);
b_eq_y(4,0) = vy_fin;
b_eq_y(5,0) = ay_fin;


b_eq_z(0,0) = X_init(2);
b_eq_z(1,0) = vz_init;
b_eq_z(2,0) = az_init;
b_eq_z(3,0) = X_fin(2);
b_eq_z(4,0) = vz_fin;
b_eq_z(5,0) = az_fin;

Eigen::MatrixXd cost_mat(17,17);
Eigen::MatrixXd cost(11,11);
Eigen::MatrixXd zero_mat( 6,6);
Eigen::MatrixXd cost_mat_inv(17,17);

Eigen::MatrixXd lincost_x(11 ,1);
Eigen::MatrixXd lincost_y(11 ,1);
Eigen::MatrixXd lincost_z(11 ,1);


Eigen::MatrixXd x_waypoints(numWayPts,1);
Eigen::MatrixXd y_waypoints(numWayPts,1);
Eigen::MatrixXd z_waypoints(numWayPts,1);

for( int l =0; l< numWayPts ; l++){
x_waypoints( l,0)= Path_vector( l , 0 );
y_waypoints( l,0)= Path_vector( l , 1 );
z_waypoints( l,0)= Path_vector( l , 2 );

}

/*
zero_mat.setZero();
float weight_smoothness = 10 ;
// std::cout << " Till Here 0 ============================" << std::endl;

cost = P_in.transpose()*P_in + weight_smoothness*Pddot_in.transpose()*Pddot_in;
// std::cout << " Till Here 1 ============================" << std::endl;

cost_mat.block( 0,0 , 11 ,11) =  cost;
// std::cout << " Till Here 1 ============================" << std::endl;

cost_mat.block( 0,11, 11 ,6) =  A_eq.transpose();
// std::cout << " Till Here 1 ============================" << std::endl;

cost_mat.block(11 ,0 , 6,11 ) = A_eq;
// std::cout << " Till Here 1 ============================" << std::endl;

cost_mat.block( 0 ,11 , 6,6) = zero_mat;


// std::cout << " Till Here 1 ============================" << std::endl;

std::cout<< P_in << std::endl;
cost_mat_inv = cost_mat.inverse();
lincost_x = -( P_in.transpose()*x_waypoints);
lincost_y = -( P_in.transpose()*y_waypoints);
lincost_z = -( P_in.transpose()*z_waypoints);

// std::cout << " Till Here 2 ============================" << std::endl;


Eigen::MatrixXf temp_matx( 17,1);
Eigen::MatrixXf temp_maty( 17,1);
Eigen::MatrixXf temp_matz( 17,1);

// std::cout << " Till Here 3 ============================" << std::endl;


temp_matx.block(0 ,0 , 11,1)=  -lincost_x;
temp_maty.block(0 ,0 , 11,1)=  -lincost_y;
temp_matz.block(0 ,0 , 11,1)=  -lincost_z;

// std::cout << " Till Here 4 ============================" << std::endl;


temp_matx.block(11 ,0 , 6,1) = b_eq_x;
temp_maty.block(11 ,0 , 6,1) = b_eq_y;
temp_matz.block(11 ,0 , 6,1) = b_eq_z;

// std::cout << " Till Here 5 ============================" << std::endl;

Eigen::MatrixXf solx( 17,1); 
Eigen::MatrixXf soly( 17,1); 
Eigen::MatrixXf solz( 17,1); 

// std::cout << " Till Here 6 ============================" << std::endl;

solx = cost_mat_inv*temp_matx;
soly = cost_mat_inv*temp_maty;
solz = cost_mat_inv*temp_matz;




Eigen::MatrixXf cx( 11, 1);
Eigen::MatrixXf cy( 11, 1);
Eigen::MatrixXf cz( 11, 1);
// std::cout << " Till Here 8 ============================" << std::endl;

cx=  solx.block( 0 , 0 , 11,1);
cy=  soly.block( 0 , 0 , 11,1);
cz=  solz.block( 0 , 0 , 11,1);

// std::cout << " Till Here 9 ============================" << std::endl;

std::vector<Eigen::MatrixXf> coefficients_vector;
coefficients_vector.push_back(cx);
coefficients_vector.push_back(cy);
coefficients_vector.push_back(cz);
*/


}
/*
std::vector<Eigen::MatrixXf> Bernstein::Bernstein_polynomial_function::get_10th_order_bernstein_matrix( float time_to_desination , int num_waypts)
{

std::cout << time_to_desination << " *******************" << std::endl;

Eigen::MatrixXf P(num_waypts , 11 ); // Bernstein matrix of 10 the order
Eigen::MatrixXf Pdot(num_waypts , 11 ); // Bernstein matrix of 10 the order
Eigen::MatrixXf Pddot(num_waypts , 11 ); // Bernstein matrix of 10 the order


Eigen::MatrixXf temp_P(num_waypts , 11 ); 
Eigen::MatrixXf temp_Pdot(num_waypts , 11 ); 
Eigen::MatrixXf temp_Pddot(num_waypts , 11 );


Eigen::MatrixXf delta_time_vector(num_waypts , 1);
std::vector<Eigen::MatrixXf> Bernstein_matrices;

float delta_t = 0.5;

for( int i=0 ; i < num_waypts ; i++){
	delta_time_vector(i , 0 ) = (time_to_desination/num_waypts )*i;
}



temp_P = get_P_matrix_10th_order(P , delta_time_vector , num_waypts);
temp_Pdot = get_Pdot_matrix_10th_order( Pdot, delta_time_vector , num_waypts);
temp_Pddot = get_Pddot_matrix_10th_order( Pddot, delta_time_vector , num_waypts);


Bernstein_matrices.push_back(P);
Bernstein_matrices.push_back(Pdot);
Bernstein_matrices.push_back(Pddot);

return Bernstein_matrices;

}
*/


void Bernstein::Bernstein_polynomial_function::get_10th_order_bernstein_matrix(float execTime , int numWayPts)
{
    float tMax = execTime;
    float tMin = 0.0;
    float l = (tMax - tMin);
    float delT = tMax/numWayPts;
    int   n = 10;

    float t = 0;

    Eigen::MatrixXd P_(numWayPts, order + 1);
    Eigen::MatrixXd Pdot_(numWayPts, order + 1);
    Eigen::MatrixXd Pddot_(numWayPts, order + 1);
    
    for(int i = 0; i<numWayPts; i++)
    {
        t += float(i)*delT;

        for(int r = 0; r<=n; r++)
        {
            float t_ = (t-tMin)/l;

            P_(i,r)  =  binomialCoefficient(n,r)*(pow((1-t_), (n-r))*pow(t_, r));

            switch(r)
            {
                case 0:
                    Pdot_(i,r)  = -10.0*pow((-t_ + 1),9);
                    Pddot_(i,r) =  90.0*pow((-t_ + 1),8);

                case 1:
                    Pdot_(i,r)  = -90.0*t_*pow((-t_ + 1),8) + 10.0*pow((-t_ + 1),9);
                    Pddot_(i,r) = 720.0*t_*pow((-t_ + 1), 7) - 180.0*pow((-t_ + 1),8);
                
                case 2:
                    Pdot_(i,r)  = -360.0*pow(t_,2)*pow((-t_ + 1),7) + 90.0*t_*pow((-t_ + 1),8);
                    Pddot_(i,r) = 2520.0*pow(t_,2)*pow((-t_ + 1),6) - 1440.0*t*pow((-t_ + 1),7) + 90.0*pow((-t_ + 1),8);
                
                case 3:
                    Pdot_(i,r)  = -840.0*pow(t_,3)*pow((-t_ + 1),6) + 360.0*pow(t_,2)*pow((-t_ + 1),7);
                    Pddot_(i,r) = 5040.0*pow(t_,3)*pow((-t_ + 1),5) - 5040.0*pow(t_,2)*pow((-t_ + 1),6) + 720.0*t_*pow((-t_ + 1),7);
                
                case 4:
                    Pdot_(i,r)  = -1260.0*pow(t_,4)*pow((-t_ + 1),5) + 840.0*pow(t_,3)*pow((-t_ + 1),6);
                    Pddot_(i,r) = 6300.0*pow(t_,4)*pow((-t_ + 1),4) - 10080.0*pow(t_,3)*pow((-t_ + 1),5) + 2520.0*pow(t_,2)*pow((-t_ + 1),6);
                
                case 5:
                    Pdot_(i,r)  = -1260.0*pow(t_,5)*pow((-t_ + 1),4) + 1260.0*pow(t_,4)*pow((-t_ + 1),5);
                    Pddot_(i,r) = 5040.0*pow(t_,5)*pow((-t_ + 1),3) - 12600.0*pow(t_,4)*pow((-t_ + 1),4) + 5040.0*pow(t_,3)*pow((-t_ + 1),5);

                case 6:
                    Pdot_(i,r)  = -840.0*pow(t_,6)*pow((-t_ + 1),3) + 1260.0*pow(t_,5)*pow((-t_ + 1),4);
                    Pddot_(i,r) = 2520.0*pow(t_,6)*pow((-t_ + 1),2) - 10080.0*pow(t_,5)*pow((-t_ + 1),3) + 6300.0*pow(t_,4)*pow((-t_ + 1),4);
                
                case 7:
                    Pdot_(i,r)  = -360.0*pow(t_,7)*pow((-t_ + 1),2) + 840.0*pow(t_,6)*pow((-t_ + 1),3);
                    Pddot_(i,r) = -360.0*pow(t_,7)*(2*t - 2) - 5040.0*pow(t_,6)*pow((-t_ + 1),2) + 5040.0*pow(t,5)*pow((-t_ + 1),3);

                case 8:
                    Pdot_(i,r)  = 45.0*pow(t_,8)*(2*t_ - 2) + 360.0*pow(t_,7)*pow((-t_ + 1),2);
                    Pddot_(i,r) = 90.0*pow(t_,8) + 720.0*pow(t_,7)*(2*t_ - 2) + 2520.0*pow(t_,6)*pow((-t_ + 1),2);

                case 9:
                    Pdot_(i,r)  = -10.0*pow(t_,9) + 9*pow(t_,8)*(-10.0*t_ + 10.0);
                    Pddot_(i,r) = -180.0*pow(t_,8) + 72*pow(t_,7)*(-10.0*t_ + 10.0);

                case 10:
                    Pdot_(i,r)  = 90.0*pow((-t_ + 1), 8);
                    Pddot_(i,r) = 90.0*pow(t_, 8);   

            }

        }
    }



    P     = &P_;
    Pdot  = &Pdot_;
    Pddot = &Pddot_;


    std::cout<<"Calculated Coefficient Matrices " <<std::endl;
}


/*
Eigen::MatrixXf Bernstein::Bernstein_polynomial_function::get_P_matrix_10th_order(  Eigen::MatrixXf  P,   Eigen::MatrixXf delta_time_vector , int num_waypts)
{


for( int i=0 ; i< num_waypts ; i++ )
{

for( int j =0 ; j < 10 ; j++){

P( i ,j ) = binomialCoefficient( 10 , j)*pow( ( 1- delta_time_vector(i ,0)) , ( 10-j) )* pow( delta_time_vector(i ,0) , j );  //((1-t)**(n-0))*t**0

}

}
return P;

}



Eigen::MatrixXf Bernstein::Bernstein_polynomial_function::get_Pdot_matrix_10th_order(  Eigen::MatrixXf  Pdot,   Eigen::MatrixXf delta_time_vector , int num_waypts)
{

for( int i=0 ; i< num_waypts ; i++ )
{
Pdot( i , 0) =  -10.0*pow( (-delta_time_vector(i ,0) + 1) , 9) ;  
Pdot( i , 1) =  -90.0*delta_time_vector(i ,0)*pow( (-delta_time_vector(i ,0) + 1) , 8) + pow( 10*(-delta_time_vector(i ,0) + 1),9) ;  //  -90.0*t*(-t + 1)**8 + 10.0*(-t + 1)**9
Pdot( i , 2) =  -360*( pow( delta_time_vector(i ,0) ,2 ))*( pow( -delta_time_vector(i ,0) ,7 )) ;  //-360.0*t**2*(-t + 1)**7 + 90.0*t*(-t + 1)**8
Pdot( i , 3) =  -840.0*pow( delta_time_vector(i ,0) , 3)*pow((-delta_time_vector(i ,0) +1) , 6) +360*pow( delta_time_vector(i ,0) , 2)*pow((-delta_time_vector(i ,0) +1) , 7) ;  
//-840.0*t**3*(-t + 1)**6 + 360.0*t**2*(-t + 1)**7
Pdot( i , 4) =  -1260.0*pow( delta_time_vector(i ,0) , 4)*pow((-delta_time_vector(i ,0) +1) , 5) + 840*pow( delta_time_vector(i ,0) , 3)*pow((-delta_time_vector(i ,0) +1) , 6) ; 
 //-1260.0*t**4*(-t + 1)**5 + 840.0*t**3*(-t + 1)**6

Pdot( i , 5) =  -1260.0*pow( delta_time_vector(i ,0) , 5)*pow((-delta_time_vector(i ,0) +1) , 4) +1260*pow( delta_time_vector(i ,0) , 4)*pow((-delta_time_vector(i ,0) +1) , 5) ;  //  -1260.0*t**5*(-t + 1)**4 + 1260.0*t**4*(-t + 1)**5

Pdot( i , 6) =  -840.0*pow( delta_time_vector(i ,0) , 6)*pow((-delta_time_vector(i ,0) +1) , 3) +1260*pow( delta_time_vector(i ,0) , 5)*pow((-delta_time_vector(i ,0) +1) , 4) ;  //  -840.0*t**6*(-t + 1)**3 + 1260.0*t**5*(-t + 1)**4

Pdot( i , 7) =  -360.0*pow( delta_time_vector(i ,0) , 7)*pow( (-delta_time_vector(i ,0) +1)  , 2) + 840*pow(delta_time_vector(i ,0) , 6)*pow( (-delta_time_vector(i ,0) +1)  , 3) ;  // -360.0*t**7*(-t + 1)**2 + 840.0*t**6*(-t + 1)**3

Pdot( i , 8) =  45.0*pow( delta_time_vector(i ,0)  , 8)*( 2*delta_time_vector(i ,0) -2 ) + 360*pow( delta_time_vector(i ,0)  , 7)*pow( (-delta_time_vector(i ,0) +1)  , 2) ;  // 45.0*t**8*(2*t - 2) + 360.0*t**7*(-t + 1)**2

Pdot( i , 9) =  -10.0*pow( delta_time_vector(i ,0) , 9) + 9*(pow( delta_time_vector(i ,0) , 8))*( -10*delta_time_vector(i ,0) +10 ) ;  //-10.0*t**9 + 9*t**8*(-10.0*t + 10.0)

Pdot( i , 10) =  10.0*pow( delta_time_vector(i ,0) , 9) ;  //10.0*t**9


}


return Pdot;

}



Eigen::MatrixXf Bernstein::Bernstein_polynomial_function::get_Pddot_matrix_10th_order(  Eigen::MatrixXf  Pddot,   Eigen::MatrixXf delta_time_vector , int num_waypts)
{
float t; 
for( int i=0 ; i< num_waypts ; i++ )
{

t = delta_time_vector(i ,0);
Pddot( i , 0) =  -90.0*pow( (-t + 1) , 8) ;  // 90.0*(-t + 1)**8
Pddot( i , 1) =  720.0*t*pow( (-t + 1) , 7) - 180*pow( (-t + 1),8) ; 
 //720.0*t*(-t + 1)**7 - 180.0*(-t + 1)**8
Pddot( i , 2) =  -2520*( pow( t,2 ))*pow( (-t+1), 6 ) -1440*t*pow (( -t+1) ,7) + 90*pow( ( -t+1) , 8 ) ; 
 //2520.0*t**2*(-t + 1)**6 - 1440.0*t*(-t + 1)**7 + 90.0*(-t + 1)**8
Pddot( i , 3) = 5040*( pow( t,3 ))*pow( (-t+1), 5 ) - 5040*pow(t ,2)*pow (( -t+1) ,6) + 720*t*pow( ( -t+1) , 7 ) ; 
//5040.0*t**3*(-t + 1)**5 - 5040.0*t**2*(-t + 1)**6 + 720.0*t*(-t + 1)**7
Pddot( i , 4) =  6300*pow( t,4)*pow((-t+1) , 4) - 10080*pow(t,3)*pow((-t+1) , 5) + 2520*pow(t,2)*pow((-t+1) , 6);  
//6300.0*t**4*(-t + 1)**4 - 10080.0*t**3*(-t + 1)**5 + 2520.0*t**2*(-t + 1)**6
Pddot( i , 5) =  5040*pow( t,5)*pow((-t+1) , 3) - 12600*pow(t,4)*pow((-t+1) , 4) + 5040*pow(t,3)*pow((-t+1) , 5) ;  
//5040.0*t**5*(-t + 1)**3 - 12600.0*t**4*(-t + 1)**4 + 5040.0*t**3*(-t + 1)**5
Pddot( i , 6) = 2520*pow( t,6)*pow((-t+1) , 2) - 10080*pow(t,5)*pow((-t+1) , 3) + 6300*pow(t,4)*pow((-t+1) , 4) ;  
//2520.0*t**6*(-t + 1)**2 - 10080.0*t**5*(-t + 1)**3 + 6300.0*t**4*(-t + 1)**4
Pddot( i , 7) = -360*pow( t,7)*(2*t-2) - 5040*pow(t,6)*pow((-t+1) , 2) + 5040*pow(t,5)*pow((-t+1) , 3) ; 
// -360.0*t**7*(2*t - 2) - 5040.0*t**6*(-t + 1)**2 + 5040.0*t**5*(-t + 1)**3
Pddot( i , 8) =  90*pow( t,8) + 720*pow(t,7)*(2*t-2)  + 2520*pow(t,6)*pow((-t+1) , 2) ;  
// 90.0*t**8 + 720.0*t**7*(2*t - 2) + 2520.0*t**6*(-t + 1)**2
Pddot( i , 9) =  -180*pow( t,8) + 72*pow(t,7)*(10*t +10) ;  
// -180.0*t**8 + 72*t**7*(-10.0*t + 10.0)
Pddot( i , 10) =  90.0*pow( t  , 8) ; 
// 90.0*t**8
}


return Pddot;

}



*/
float Bernstein::Bernstein_polynomial_function::binomialCoefficient(int n, int i)
{
    return float(factorial(n))/(float(factorial(n-i))*float(factorial(i)));
}



float Bernstein::Bernstein_polynomial_function::factorial(int x)
{
    if(x==1 || x==0)
    {
        return 1;
    }

    else
    {
        return float(x*factorial(x-1));
    }
}