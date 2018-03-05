#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = (3 - n_aug_);

  weights_ = VectorXd(2*n_aug_+1);


  // Generation of sigma points
  x_aug_ = VectorXd(n_aug_);
  P_aug_ = MatrixXd(n_aug_, n_aug_);
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);


  // Prediction of sigma points
  Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);



  // Measurement space dimension (Laser)
  n_z_radar = 3;
  // Measurement space dimension (Laser)
  n_z_laser = 2;

  // Transform sigma points into measurement space (Radar)
  Zsig_radar = MatrixXd(n_z_radar, 2 * n_aug + 1);

  // Transform sigma points into measurement space (Laser)
  Zsig_laser = MatrixXd(n_z_laser, 2 * n_aug + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  // set weights
  double weight_0 = lambda/(lambda+n_aug);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug+lambda);
    weights_(i) = weight;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  // Prediction

  // Step-1 : Generate Sigma points Matrix

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //Update augmented sigma points
  Xsig_aug_.col(0)  = x_aug_;

  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug_.col(i + 1)          = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
  }


  // Step-2 : Use process model function to Predict sigma points
  for (int i = 0; i < (2 * n_aug + 1); i++)
  {
    
    double p_x  = Xsig_aug(0, i);
    double p_y  = Xsig_aug(1, i);
    double v    = Xsig_aug(2, i);
    double yaw  = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // predicted state values
    double px_pred, py_pred;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_pred = p_x + v/yawd * ( sin(yaw + (yawd * delta_t)) - sin(yaw));
        py_pred = p_y + v/yawd * ( cos(yaw) - cos(yaw + (yawd * delta_t)) );
    }
    else {
        px_pred = p_x + (v * delta_t * cos(yaw));
        py_pred = p_y + (v * delta_t * sin(yaw));
    }

    double v_pred = v;
    double yaw_pred = yaw + (yawd * delta_t);
    double yawd_pred = yawd;

    // add noise values
    px_pred = px_pred + (0.5 * nu_a * delta_t * delta_t * cos(yaw));
    py_pred = py_pred + (0.5 * nu_a * delta_t * delta_t * sin(yaw));
    v_pred = v_pred + (nu_a * delta_t);

    yaw_pred = yaw_pred + (0.5 * nu_yawdd * delta_t * delta_t);
    yawd_pred = yawd_pred + (nu_yawdd * delta_t);

    //write predicted sigma point into right column
    Xsig_pred(0, i) = px_pred;
    Xsig_pred(1, i) = py_pred;
    Xsig_pred(2, i) = v_pred;
    Xsig_pred(3, i) = yaw_pred;
    Xsig_pred(4, i) = yawd_pred;
  }


  // Step - 3: Predict the mean and covariance matrix

  //predicted state mean
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {
    x_ = x_ + weights_(i) * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  for (int i = 0; i < 2 * n_aug + 1; i++) {

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;

    //angle normalization
    while (x_diff(3) > M_PI) {
      x_diff(3) -= (2. * M_PI);
    }

    while (x_diff(3) < -M_PI) {
      x_diff(3) += (2. * M_PI);
    }

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}



/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}



/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // Step - 1 : Transform the predicted mean and state covirance (5x1 matrix) into measurement space (3x1 matrix)
  // The measurement space mean and measurement covariance matrix
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {

    double p_x = Xsig_pred(0, i);
    double p_y = Xsig_pred(1, i);
    double v   = Xsig_pred(2, i);
    double yaw = Xsig_pred(3, i);

    double v1 = (cos(yaw) * v);
    double v2 = (sin(yaw) * v);

    // measurement model
    //rho
    Zsig_radar(0, i) = sqrt((p_x * p_x) + (p_y * p_y));
    //phi
    Zsig_radar(1, i) = atan2(p_y, p_x);
    //rho_dot
    Zsig_radar(2, i) = ((p_x * v1) + (p_y * v2) ) / sqrt((p_x * p_x) + (p_y * p_y));
  }


  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  for (int i=0; i < (2 * n_aug + 1); i++) {
      z_pred = z_pred + weights_(i) * Zsig_radar.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  for (int i = 0; i < (2 * n_aug + 1); i++) {

    //residual
    VectorXd z_diff = Zsig_radar.col(i) - z_pred;

    //angle normalization
    while (z_diff(1) > M_PI) {
      z_diff(1) -= (2. * M_PI);
    }

    while (z_diff(1) < -M_PI) {
      z_diff(1) += (2. * M_PI);
    }

    S = S + (weights_(i) * z_diff * z_diff.transpose());

  }


  // Add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R <<  (std_radr_ * std_radr_), 0, 0,
        0, (std_radphi_ * std_radphi_), 0,
        0, 0, (std_radrd_ * std_radrd_);


  S = S + R;




  // Step - 2 : Use the actual measurement values to update the state mean and covariance matrix
  

}
