#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse = VectorXd(4);
  rmse << 0, 0, 0, 0;

  // If estimations size is not equal to ground truth, cannot calculate rmse
  if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
  	cout<<"Invalid estimation or ground truth data size"<<endl;
  	return rmse;
  }

  for(unsigned int i = 0; i < estimations.size(); ++i) {
  	// Get the residual (Error in RMSE)
  	VectorXd residual = (estimations[i] - ground_truth[i]);

  	// Get the square of residual (Squared in RMSE)
  	residual = (residual.array() * residual.array());

  	//Add the residual to overall rmse (Add the squared errors)
  	rmse += residual;
  }

  // Get the average of all rmse values (Mean in RMSE)
  rmse = rmse/estimations.size();

  //Calculate the square root of summed rmse (Root mean squared error)
  rmse = rmse.array().sqrt();

  return rmse;
}