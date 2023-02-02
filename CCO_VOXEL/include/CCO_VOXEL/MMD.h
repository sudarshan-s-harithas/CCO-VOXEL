
#include "utils.h"
#include <random>
#include <algorithm>
#include <fstream>
namespace MMDFunctions
{

  class MMD_variants
  {

  public:
    Eigen::Matrix<float, 100, 5> Weights;

    float MMD_vectorized(Eigen::MatrixXf actual_distribution);
    float MMD_interpolation_method(float dist);
    std::string assign_weights(std::string path_to_weights);
    float MMD_transformed_features(Eigen::MatrixXf actual_distribution);
    float MMD_transformed_features_RBF(Eigen::MatrixXf actual_distribution);
    float RBF_kernel(float val1, float val2);
  };

}

std::string MMDFunctions::MMD_variants::assign_weights(std::string path_to_weights)
{
  int rows = 100;
  int cols = 5;
  std::string file = path_to_weights;
  std::ifstream in(file);

  std::string line;

  int row = 0;
  int col = 0;

  Eigen::MatrixXf res = Eigen::MatrixXf(rows, cols);

  if (in.is_open())
  {

    while (std::getline(in, line))
    {

      char *ptr = (char *)line.c_str();
      int len = line.length();

      col = 0;

      char *start = ptr;
      for (int i = 0; i < len; i++)
      {

        if (ptr[i] == ',')
        {
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

float MMDFunctions::MMD_variants::MMD_vectorized(Eigen::MatrixXf actual_distribution)
{

  int num_samples_of_distance_distribution = actual_distribution.cols();
  Eigen::MatrixXf squared_dist(1, num_samples_of_distance_distribution);
  Eigen::MatrixXf temp(1, num_samples_of_distance_distribution);
  Eigen::MatrixXf alpha_weights(1, num_samples_of_distance_distribution);

  alpha_weights.setOnes();

  // temp = actual_distribution.array();

  squared_dist = (actual_distribution.array()).square();

  Eigen::MatrixXf kernel_matrix1(num_samples_of_distance_distribution, num_samples_of_distance_distribution);
  Eigen::MatrixXf temp_kernel_matrix1(num_samples_of_distance_distribution, num_samples_of_distance_distribution);
  Eigen::MatrixXf temp_kernel_matrix2(num_samples_of_distance_distribution, num_samples_of_distance_distribution);

  Eigen::MatrixXf One_matrix(num_samples_of_distance_distribution, num_samples_of_distance_distribution);
  One_matrix.setOnes();

  // std::cout << " Entered the main funciton /////////////" << std::endl;
  temp_kernel_matrix1 = actual_distribution.transpose() * actual_distribution;

  temp_kernel_matrix2 = squared_dist.transpose() * squared_dist;

  kernel_matrix1 = One_matrix + 2 * temp_kernel_matrix1 + temp_kernel_matrix2;

  Eigen::Matrix<float, 1, 1> res;
  Eigen::Matrix<float, 1, 1> res1;
  Eigen::Matrix<float, 1, 1> res2;
  Eigen::Matrix<float, 1, 1> res3;

  float temp_value = 0;

  double res_val;
  res1 = alpha_weights * kernel_matrix1 * (alpha_weights.transpose());

  res2(0, 0) = num_samples_of_distance_distribution;
  res3(0, 0) = num_samples_of_distance_distribution;

  res(0, 0) = res1(0, 0) - 2 * res2(0, 0) + res3(0, 0);
  res_val = double(res(0, 0));

  return float(res_val);
}

float MMDFunctions::MMD_variants::MMD_interpolation_method(float dist)
{

  Eigen::MatrixXf coefficient_matrix(1, 7);
  Eigen::MatrixXf distance_matrix(1, 7);
  Eigen::MatrixXf result(1, 1);

  coefficient_matrix(0, 0) = 7720.61614034;
  coefficient_matrix(0, 1) = -12367.93260857;
  coefficient_matrix(0, 2) = -1541.47268526;
  coefficient_matrix(0, 3) = 14365.74417541;
  coefficient_matrix(0, 4) = -10918.17820688;
  coefficient_matrix(0, 5) = 3339.69565025;
  coefficient_matrix(0, 6) = -373.83984028;

  for (int i = 0; i < 7; i++)
  {
    distance_matrix(0, i) = pow(dist, i);
  }

  result = coefficient_matrix * distance_matrix.transpose();

  return result(0, 0);
}

float MMDFunctions::MMD_variants::MMD_transformed_features(Eigen::MatrixXf actual_distribution)
{

  Eigen::MatrixXf transformed_features(1, 5);

  transformed_features = actual_distribution * Weights;

  // std::cout << transformed_features << std::endl;

  int num_samples_of_distance_distribution = transformed_features.cols();
  Eigen::MatrixXf squared_dist(1, num_samples_of_distance_distribution);
  Eigen::MatrixXf temp(1, num_samples_of_distance_distribution);
  Eigen::MatrixXf alpha_weights(1, num_samples_of_distance_distribution);

  alpha_weights.setOnes();

  // temp = actual_distribution.array();

  squared_dist = (transformed_features.array()).square();

  Eigen::MatrixXf kernel_matrix1(num_samples_of_distance_distribution, num_samples_of_distance_distribution);
  Eigen::MatrixXf temp_kernel_matrix1(num_samples_of_distance_distribution, num_samples_of_distance_distribution);
  Eigen::MatrixXf temp_kernel_matrix2(num_samples_of_distance_distribution, num_samples_of_distance_distribution);

  Eigen::MatrixXf One_matrix(num_samples_of_distance_distribution, num_samples_of_distance_distribution);
  One_matrix.setOnes();

  // std::cout << " Entered the main funciton /////////////" << std::endl;
  temp_kernel_matrix1 = transformed_features.transpose() * transformed_features;

  temp_kernel_matrix2 = squared_dist.transpose() * squared_dist;

  kernel_matrix1 = One_matrix + 2 * temp_kernel_matrix1 + temp_kernel_matrix2;

  Eigen::Matrix<float, 1, 1> res;
  Eigen::Matrix<float, 1, 1> res1;
  Eigen::Matrix<float, 1, 1> res2;
  Eigen::Matrix<float, 1, 1> res3;

  float temp_value = 0;

  double res_val;
  res1 = alpha_weights * kernel_matrix1 * (alpha_weights.transpose());

  res2(0, 0) = num_samples_of_distance_distribution;
  res3(0, 0) = num_samples_of_distance_distribution;

  res(0, 0) = res1(0, 0) - 2 * res2(0, 0) + res3(0, 0);
  res_val = double(res(0, 0));

  return float(res_val);
}

float MMDFunctions::MMD_variants::RBF_kernel(float val1, float val2)
{

  float data1 = std::pow((val1 - val2), 2) / (-2 * std::pow(0.1, 2));
  float data2 = std::exp(data1);

  return data2;
}

float MMDFunctions::MMD_variants::MMD_transformed_features_RBF(Eigen::MatrixXf actual_distribution)
{

  Eigen::MatrixXf transformed_features(1, 5);

  transformed_features = actual_distribution * Weights;

  // std::cout << transformed_features << std::endl;

  int num_samples_of_distance_distribution = transformed_features.cols();

  Eigen::MatrixXf Matrix1(num_samples_of_distance_distribution, num_samples_of_distance_distribution);
  Eigen::MatrixXf Matrix2(num_samples_of_distance_distribution, num_samples_of_distance_distribution);
  Eigen::MatrixXf Matrix3(num_samples_of_distance_distribution, num_samples_of_distance_distribution);

  for (int i = 0; i < num_samples_of_distance_distribution; i++)
  {

    for (int j = 0; j < num_samples_of_distance_distribution; j++)
    {

      Matrix1(i, j) = RBF_kernel(transformed_features(0, i), transformed_features(0, j));
      Matrix2(i, j) = RBF_kernel(transformed_features(0, i), 0.0);
      Matrix3(i, j) = RBF_kernel(0.0, 0.0);
    }
  }

  Eigen::MatrixXf alpha_weights(1, num_samples_of_distance_distribution);
  alpha_weights.setOnes();

  Eigen::Matrix<float, 1, 1> cost1;
  Eigen::Matrix<float, 1, 1> cost2;
  Eigen::Matrix<float, 1, 1> cost3;

  cost1 = alpha_weights * Matrix1 * alpha_weights.transpose();
  cost2 = alpha_weights * Matrix2 * alpha_weights.transpose();
  cost3 = alpha_weights * Matrix3 * alpha_weights.transpose();

  float MMD_RBF_cost = cost1(0, 0) - 2 * cost2(0, 0) + cost3(0, 0);

  return MMD_RBF_cost;
}
