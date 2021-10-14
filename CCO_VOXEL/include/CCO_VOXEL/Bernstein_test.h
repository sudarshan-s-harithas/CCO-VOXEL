#include <random>
#include"utils.h"

namespace Bernstein {

class Bernstein_polynomial_function{

public:

	int order =10;
	std::vector<Eigen::Vector3d> coeffs;

	std::vector<Eigen::Vector3d> get_coefficients(std::vector<Eigen::Vector3d> wayPts  , int numWayPts , float execTime );
    std::vector<Eigen::Vector3d> perturb_coefficients( std::vector<Eigen::Vector3d> Bernstein_coeff );
	float binomial(int n, int i);
	float factorial(int x);





};

}




std::vector<Eigen::Vector3d> Bernstein::Bernstein_polynomial_function::perturb_coefficients( std::vector<Eigen::Vector3d> Bernstein_coeff )
{



    
}


std::vector<Eigen::Vector3d>  Bernstein::Bernstein_polynomial_function::get_coefficients( std::vector<Eigen::Vector3d> wayPts  , int numWayPts ,  float execTime )
{



    float tMax = execTime;
    float tMin = 0.0;
    float l = (tMax - tMin);
    float delT = tMax/numWayPts;
    int   n = order;

    float t = 0;

    Eigen::MatrixXd P(numWayPts, order + 1);
    Eigen::MatrixXd Pdot(numWayPts, order + 1);
    Eigen::MatrixXd Pddot(numWayPts, order + 1);
    for(int i = 0; i<numWayPts; i++)
    {
        float t_ = float(float(i)/float(numWayPts-1));
        //t_  = t_/l;

        std::cout<<"t_ is: "<<t_<<std::endl;

        for(int r = 0; r<=n; r++)
        {
            P(i,r)  =  binomial(n,r) * pow((1-t_), (n-r)) * pow(t_, r);

            switch(r)
            {
                case 0:
                    Pdot(i,r)  = -10.0*pow((-t_ + 1),9);
                    Pddot(i,r) =  90.0*pow((-t_ + 1),8);
                    break;

                case 1:
                    Pdot(i,r)  = -90.0*t_*pow((-t_ + 1),8) + 10.0*pow((-t_ + 1),9);
                    Pddot(i,r) = 720.0*t_*pow((-t_ + 1), 7) - 180.0*pow((-t_ + 1),8);
                    break;
                
                case 2:
                    Pdot(i,r)  = -360.0*pow(t_,2)*pow((-t_ + 1),7) + 90.0*t_*pow((-t_ + 1),8);
                    Pddot(i,r) = 2520.0*pow(t_,2)*pow((-t_ + 1),6) - 1440.0*t_*pow((-t_ + 1),7) + 90.0*pow((-t_ + 1),8);
                    break;
                
                case 3:
                    Pdot(i,r)  = -840.0*pow(t_,3)*pow((-t_ + 1),6) + 360.0*pow(t_,2)*pow((-t_ + 1),7);
                    Pddot(i,r) = 5040.0*pow(t_,3)*pow((-t_ + 1),5) - 5040.0*pow(t_,2)*pow((-t_ + 1),6) + 720.0*t_*pow((-t_ + 1),7);
                    break;
                
                case 4:
                    Pdot(i,r)  = -1260.0*pow(t_,4)*pow((-t_ + 1),5) + 840.0*pow(t_,3)*pow((-t_ + 1),6);
                    Pddot(i,r) = 6300.0*pow(t_,4)*pow((-t_ + 1),4) - 10080.0*pow(t_,3)*pow((-t_ + 1),5) + 2520.0*pow(t_,2)*pow((-t_ + 1),6);
                    break;
                
                case 5:
                    Pdot(i,r)  = -1260.0*pow(t_,5)*pow((-t_ + 1),4) + 1260.0*pow(t_,4)*pow((-t_ + 1),5);
                    Pddot(i,r) = 5040.0*pow(t_,5)*pow((-t_ + 1),3) - 12600.0*pow(t_,4)*pow((-t_ + 1),4) + 5040.0*pow(t_,3)*pow((-t_ + 1),5);
                    break;

                case 6:
                    Pdot(i,r)  = -840.0*pow(t_,6)*pow((-t_ + 1),3) + 1260.0*pow(t_,5)*pow((-t_ + 1),4);
                    Pddot(i,r) = 2520.0*pow(t_,6)*pow((-t_ + 1),2) - 10080.0*pow(t_,5)*pow((-t_ + 1),3) + 6300.0*pow(t_,4)*pow((-t_ + 1),4);
                    break;
                
                case 7:
                    Pdot(i,r)  = -360.0*pow(t_,7)*pow((-t_ + 1),2) + 840.0*pow(t_,6)*pow((-t_ + 1),3);
                    Pddot(i,r) = -360.0*pow(t_,7)*(2*t_ - 2) - 5040.0*pow(t_,6)*pow((-t_ + 1),2) + 5040.0*pow(t_,5)*pow((-t_ + 1),3);
                    break;

                case 8:
                    Pdot(i,r)  = 45.0*pow(t_,8)*(2*t_ - 2) + 360.0*pow(t_,7)*pow((-t_ + 1),2);
                    Pddot(i,r) = 90.0*pow(t_,8) + 720.0*pow(t_,7)*(2*t_ - 2) + 2520.0*pow(t_,6)*pow((-t_ + 1),2);
                    break;

                case 9:
                    Pdot(i,r)  = -10.0*pow(t_,9) + 9*pow(t_,8)*(-10.0*t_ + 10.0);
                    Pddot(i,r) = -180.0*pow(t_,8) + 72*pow(t_,7)*(-10.0*t_ + 10.0);
                    break;

                case 10:
                    Pdot(i,r)  = 10.0*pow(t_, 9);//90.0*pow((-t_ + 1), 8);
                    Pddot(i,r) = 90.0*pow(t_, 8);   
                    break;


            }

    

                Pdot(i,r) = Pdot(i,r)/l;
                Pddot(i,r) = Pddot(i,r)/pow(l,2);

        }

        t = t_;

    }



Eigen::MatrixXd P_in(numWayPts , 11);
 Eigen::MatrixXd Pdot_in(numWayPts,11);
 Eigen::MatrixXd Pddot_in(numWayPts,11);

P_in =P;
Pdot_in = Pdot;
Pddot_in = Pddot;


  int numPts = wayPts.size();

    std::cout<<P<<std::endl;

    Eigen::MatrixXd xPts(numPts,1);
    Eigen::MatrixXd yPts(numPts,1);
    Eigen::MatrixXd zPts(numPts,1);


    for(int i = 0; i<wayPts.size(); i++)
    {
        Eigen::Vector3d pt = wayPts.at(i);

        xPts(i,0) = pt(0);
        yPts(i,0) = pt(1);
        zPts(i,0) = pt(2);
    }

    Eigen::Vector3d initPt = wayPts.at(0);
    Eigen::Vector3d endPt  = wayPts.at(numPts-1);

    Eigen::Vector3d initVel(0.0,0.0,0.0);
    Eigen::Vector3d endVel(0.0,0.0,0.0);
    Eigen::Vector3d initAcc(0.0,0.0,0.0);
    Eigen::Vector3d endAcc(0.0,0.0,0.0);

    // generate the matrix A_eq, B_eq_x, B_eq_y, B_eq_z
    Eigen::MatrixXd pInit(1, order+1), pdotInit(1, order+1), pddotInit(1, order+1);
    Eigen::MatrixXd pFin(1, order+1), pdotFin(1, order+1), pddotFin(1, order+1);

    int nRows = P.rows() - 1;


    for(int i = 0; i<=order; i++)
    {
        pInit(0,i) = P(0,i);
        pFin(0,i)  = P(nRows, i);

        pdotInit(0,i) = Pdot(0, i);
        pdotFin(0,i)  = Pdot(nRows, i);

        pddotInit(0, i) = Pddot(0, i);
        pddotFin(0, i)  = Pddot(nRows, i);
    }


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



b_eq_x(0,0) = initPt(0);
b_eq_x(1,0) = vx_init;
b_eq_x(2,0) = ax_init;
b_eq_x(3,0) = endPt(0);
b_eq_x(4,0) = vy_fin;
b_eq_x(5,0) = ay_fin;

b_eq_y(0,0) = initPt(1);
b_eq_y(1,0) = vy_init;
b_eq_y(2,0) = ay_init;
b_eq_y(3,0) = endPt(1);
b_eq_y(4,0) = vy_fin;
b_eq_y(5,0) = ay_fin;


b_eq_z(0,0) = initPt(2);
b_eq_z(1,0) = vz_init;
b_eq_z(2,0) = az_init;
b_eq_z(3,0) = endPt(2);
b_eq_z(4,0) = vz_fin;
b_eq_z(5,0) = az_fin;


    // std::cout<<A_eq<<std::endl;

    float weightSmoothness = 10.0;

    Eigen::MatrixXd cost = (P.transpose() * (P))  + weightSmoothness*(Pddot.transpose() * (Pddot));

    Eigen::MatrixXd costMat     = Eigen::MatrixXd::Zero(17,17);
    Eigen::MatrixXd costMatUp   = Eigen::MatrixXd::Zero(11,17);
    Eigen::MatrixXd costMatDown = Eigen::MatrixXd::Zero(6,17);

    costMatUp   << cost, A_eq.transpose();
    costMatDown << A_eq, Eigen::MatrixXd::Zero(A_eq.rows(), A_eq.rows());
    costMat     << costMatUp, costMatDown;

    // Eigen::MatrixXd costMatInv = costMat.inverse();

    Eigen::MatrixXd linCostX = -(P.transpose()*xPts);
    Eigen::MatrixXd linCostY = -(P.transpose()*yPts);
    Eigen::MatrixXd linCostZ = -(P.transpose()*zPts);


    Eigen::MatrixXd linBx(linCostX.rows() + b_eq_x.rows(),1);
    linBx << -linCostX, b_eq_x;
    
    Eigen::MatrixXd linBy(linCostY.rows() + b_eq_y.rows(),1);
    linBy << -linCostY, b_eq_x;

    Eigen::MatrixXd linBz(linCostZ.rows() + b_eq_z.rows(),1);
    linBz << -linCostZ, b_eq_x;

    Eigen::MatrixXd solX = costMat.completeOrthogonalDecomposition().solve(linBx);
    Eigen::MatrixXd solY = costMat.completeOrthogonalDecomposition().solve(linBy);
    Eigen::MatrixXd solZ = costMat.completeOrthogonalDecomposition().solve(linBz);  


    for(int i = 0; i<order+1; i++)
    {
        Eigen::Vector3d pt(solX(i,0), solY(i,0), solZ(i,0));
        std::cout<<"Coefficients are "<<pt(0)<<"\t"<<pt(1)<<"\t"<<pt(2)<<std::endl;
        coeffs.push_back(pt);
    }



return coeffs;
}



float Bernstein::Bernstein_polynomial_function::binomial(int n, int i)
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