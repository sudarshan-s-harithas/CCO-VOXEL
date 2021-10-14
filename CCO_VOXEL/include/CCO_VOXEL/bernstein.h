/** generate bernstein coefficients from fastplanner trajectory **/

#include"utils.h"
#include<random>

namespace Bernstein
{
    class BernsteinPath
    {
        public:
            int order;
            Eigen::MatrixXd P, Pdot, Pddot;             // matrix coefficients
            std::vector<Eigen::Vector3d> waypoints;
            std::vector<Eigen::Vector3d> coeffs;

            BernsteinPath();                          // default constructor
            BernsteinPath(int order_);                // parameterized constructor to change the order of beizer curve
            
            inline void generateCoeffMatrices(int numWayPts, float execTime); // generate coefficient matrices (P, Pdot, Pddot)
            inline void generateTrajCoeffs(std::vector<Eigen::Vector3d> waypts);// generate the initial coefficients using fastplanner waypoints
            inline std::vector<Eigen::MatrixXd> generatePerturbedCoeffs(int numSamples, std::vector<Eigen::Vector3d> coeffs_, Eigen::Vector3d var_vector); // generate perturbed coefficients using the current coefficients

        private:
            inline double binomial(int n, int r);
            double factorial(int n);
            inline std::vector<Eigen::MatrixXd> generateRandomCoeffs(int numSamples, std::vector<Eigen::Vector3d> coeffs_, Eigen::Vector3d var_vector); // returns vector of 3 matrices of 100x11(for x,y and z)
            // std::vector<Eigen::MatrixXd> perturb_coefficeints_of_mean_trajectory( int numSamples, std::vector<Eigen::Vector3d> mean_coeff,  Eigen::Vector3d var_vector );
    };   

}

/*--------------------------------------------- Function definitions --------------------------------------------*/


/******************************************************
* Default constructor setting order of bernstein as 0
******************************************************/
Bernstein::BernsteinPath::BernsteinPath()
{
    std::cout<<"Bernstein Initialized (default) with order as 10 "<<std::endl;
    order = 10;
}

/*****************************************************************
* Parameterized constructor setting order of bernstein as per user
******************************************************************/
Bernstein::BernsteinPath::BernsteinPath(int order_)
{
    std::cout<<"Bernstein Initialized  with order as "<<order_<<std::endl;
    order = order_;
}


/****************************************************************
* Generate bernstein coefficient matrices for given time interval
*****************************************************************/
inline void Bernstein::BernsteinPath::generateCoeffMatrices(int numWayPts, float execTime)
{
    float tMax = execTime;
    float tMin = 0.0;
    float l = (tMax - tMin);
    float delT = tMax/numWayPts;
    // std::cout<<"l and delT are "<<l<<" "<<delT<<std::endl;
    int   n = order;

    float t = 0;

    P = Eigen::MatrixXd::Zero(numWayPts, order+1);
    Pdot = Eigen::MatrixXd::Zero(numWayPts, order+1);
    Pddot = Eigen::MatrixXd::Zero(numWayPts, order + 1);
    
    for(int i = 0; i<numWayPts; i++)
    {
        float t_ = float(float(i)/float(numWayPts-1));
        //t_  = t_/l;

        // std::cout<<"t_ is: "<<t_<<std::endl;

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

}


/********************************************************************
* Generate bernstein trajectory coefficients for given time interval
*********************************************************************/
inline void Bernstein::BernsteinPath::generateTrajCoeffs(std::vector<Eigen::Vector3d> wayPts)
{
    int numPts = wayPts.size();

    // std::cout<<P<<std::endl;

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


    Eigen::MatrixXd A_eq(6,order+1);
    Eigen::MatrixXd B_eq_x(6,1), B_eq_y(6,1), B_eq_z(6,1);


    A_eq   << pInit, pdotInit, pddotInit, pFin, pdotFin, pddotFin;
    B_eq_x << initPt(0), initVel(0), initAcc(0), endPt(0), endVel(0), endAcc(0);
    B_eq_y << initPt(1), initVel(1), initAcc(1), endPt(1), endVel(1), endAcc(1);
    B_eq_z << initPt(2), initVel(2), initAcc(2), endPt(2), endVel(2), endAcc(2);

    // std::cout<<A_eq<<std::endl;

    float weightSmoothness = 30.0;

    Eigen::MatrixXd cost = (P.transpose() * (P))  + weightSmoothness*(Pddot.transpose() * (Pddot));

    Eigen::MatrixXd costMat     = Eigen::MatrixXd::Zero(17,17);
    Eigen::MatrixXd costMatUp   = Eigen::MatrixXd::Zero(11,17);
    Eigen::MatrixXd costMatDown = Eigen::MatrixXd::Zero(6,17);

    costMatUp   << cost, A_eq.transpose();
    costMatDown << A_eq, Eigen::MatrixXd::Zero(A_eq.rows(), A_eq.rows());
    costMat     << costMatUp, costMatDown;

    Eigen::MatrixXd costMatInv = costMat.inverse();

    Eigen::MatrixXd linCostX = -(P.transpose()*xPts);
    Eigen::MatrixXd linCostY = -(P.transpose()*yPts);
    Eigen::MatrixXd linCostZ = -(P.transpose()*zPts);


    Eigen::MatrixXd linBx(linCostX.rows() + B_eq_x.rows(),1);
    linBx << -linCostX, B_eq_x;
    
    Eigen::MatrixXd linBy(linCostY.rows() + B_eq_y.rows(),1);
    linBy << -linCostY, B_eq_y;

    Eigen::MatrixXd linBz(linCostZ.rows() + B_eq_z.rows(),1);
    linBz << -linCostZ, B_eq_z;
    
    // std::cout<<linBx.transpose()<<linBy.transpose()<<linBz.transpose()<<std::endl;
    // std::cout<<linBx.rows()<<" "<<linBx.cols()<<std::endl;
    Eigen::MatrixXd solX = costMatInv * linBx;
    Eigen::MatrixXd solY = costMatInv * linBy;
    Eigen::MatrixXd solZ = costMatInv * linBz; 
    coeffs.clear(); 


    for(int i = 0; i<order+1; i++)
    {
        Eigen::Vector3d pt(solX(i,0), solY(i,0), solZ(i,0));
        // std::cout<<"Coefficients are "<<pt(0)<<"\t"<<pt(1)<<"\t"<<pt(2)<<std::endl;
        coeffs.push_back(pt);
    }

    std::cout<<"Generated the coefficients"<<std::endl;

}

/*********************************************
* Function to generate random coefficients  *
**********************************************/
inline std::vector<Eigen::MatrixXd> Bernstein::BernsteinPath::generatePerturbedCoeffs(int numSamples, std::vector<Eigen::Vector3d> wayPts, Eigen::Vector3d var_vector)
{
    std::vector<Eigen::MatrixXd> perturbedCoeffs = generateRandomCoeffs(numSamples, wayPts,  var_vector);

    // now set the initial and final coordinates of each sample equal to that of initial wayPts
    return perturbedCoeffs;
}


/********************************************
* Function to calculate binomial coefficient
*********************************************/
inline double Bernstein::BernsteinPath::binomial(int n, int r)
{
    return (factorial(n)/(factorial(n-r)*factorial(r)));
}

/*********************************************
* Function to calculate factorial of a number
**********************************************/
double Bernstein::BernsteinPath::factorial(int n)
{
    if(n==1 || n == 0)
    {
        return 1.0;
    }

    if(n<0)
    {
        std::cout<<"Something wrong with factorial !!! "<<std::endl;
        return 0;
    }

    else
    {
        return double(n*factorial(n-1));
    }

}


/*******************************************
* Function to generate random trajectories *
********************************************/
inline std::vector<Eigen::MatrixXd> Bernstein::BernsteinPath::generateRandomCoeffs(int numSamples, std::vector<Eigen::Vector3d> coefficients,  Eigen::Vector3d var_vector)
{
    std::default_random_engine de(time(0));
    
    std::vector<Eigen::MatrixXd> perturbCoeffs;

    Eigen::MatrixXd coeffX(numSamples, coefficients.size()), coeffY(numSamples, coefficients.size()), coeffZ(numSamples, coefficients.size()); // 100 x 11 coefficient matrix  

    std::cout << " this size " << coefficients.size() << std::endl;
    
    for(int i = 0; i<coefficients.size(); i++)
     {
        Eigen::Vector3d pt = coefficients.at(i);

         // generate random samples around this point

        for(int j = 0; j<numSamples; j++)
        {

            // std::cout << i <<  " " << j << std::endl; 
            std::normal_distribution<double> ndX(pt(0), var_vector.x());
            std::normal_distribution<double> ndY(pt(1), var_vector.y());
            std::normal_distribution<double> ndZ(pt(2), var_vector.z());

            coeffX(j,i) = ndX(de);
            coeffY(j,i) = ndY(de);
            coeffZ(j,i) = ndZ(de);

            if(i==0 || i== 1 || i ==2 || i == coefficients.size()-1  || i == coefficients.size()-2 || i == coefficients.size()-3) 
            {
                coeffX(j,i) = pt(0);
                coeffY(j,i) = pt(1);
                coeffZ(j,i) = pt(2);
                            
            }

        }    
         
     }

     perturbCoeffs.push_back(coeffX);
     perturbCoeffs.push_back(coeffY);
     perturbCoeffs.push_back(coeffZ);

     return perturbCoeffs;
}