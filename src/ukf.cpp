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

  // Is the UKF initialized?
  is_initialized_ = false;

  time_us_ = 0;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2;
  
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

  nis_radar = 0.0;
  nis_lidar = 0.0;
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


  if (!is_initialized_) {

    // set weights
    double weight_0 = lambda_/(lambda_ + n_aug_);
    weights_(0) = weight_0;

    for (int i=1; i< (2 * n_aug_ + 1); i++) {
      double weight = 0.5/(n_aug_ + lambda_);
      weights_(i) = weight;
    }

    x_.fill(0.0);
    //x_ << meas_package.raw_measurements_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Get the raw measurements in polar coordinates
      float rho = meas_package.raw_measurements_(0);
      float theta = meas_package.raw_measurements_(1);
      float rho_dot = meas_package.raw_measurements_(2);

      x_(0) = rho * cos(theta); //px
      x_(1) = rho * sin(theta); //py
      x_(2) = rho_dot * cos(theta); //vx
      x_(3) = rho_dot * sin(theta); // vy
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // Get the raw measurements
      x_(0) = meas_package.raw_measurements_(0); //px
      x_(1) = meas_package.raw_measurements_(1); //py
    }


    P_.fill(0.0);

    for(int i = 0; i < n_x_; i++) {
      P_(i,i) = 1;
    }

    // store the starting timestamp
    time_us_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;

    std::cout<<"Done Initializing!"<<endl;
    return;
  }


  nis_radar = 0.0;
  nis_lidar = 0.0;
  
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // Delta_T in seconds
  float delta_time = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_time);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
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

  //create augmented mean state
  x_aug_.head(n_x_) = x_;
  x_aug_(5) = 0.0;
  x_aug_(6) = 0.0;

  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = (std_a_ * std_a_);
  P_aug_(6,6) = (std_yawdd_ * std_yawdd_);


  //create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();

  //Update augmented sigma points
  Xsig_aug_.col(0)  = x_aug_;

  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug_.col(i + 1)          = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
  }


  // Step-2 : Use process model function to Predict sigma points
  for (int i = 0; i < (2 * n_aug_ + 1); i++)
  {
    
    double p_x  = Xsig_aug_(0, i);
    double p_y  = Xsig_aug_(1, i);
    double v    = Xsig_aug_(2, i);
    double yaw  = Xsig_aug_(3, i);
    double yawd = Xsig_aug_(4, i);
    double nu_a = Xsig_aug_(5, i);
    double nu_yawdd = Xsig_aug_(6, i);

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
  x_.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {
    x_ = x_ + weights_(i) * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {

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


  // Transform sigma points into measurement space (Laser)
  MatrixXd Zsig_laser = MatrixXd(n_z_laser, 2 * n_aug_ + 1);


  // Step - 1 : Transform the predicted mean and state covirance (5x1 matrix) into measurement space (2x1 matrix)
  // The measurement space mean and measurement covariance matrix
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {

    double p_x = Xsig_pred(0, i);
    double p_y = Xsig_pred(1, i);

    // measurement model
    Zsig_laser(0, i) = p_x;
    Zsig_laser(1, i) = p_y;
  }


  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_laser);
  z_pred.fill(0.0);

  for (int i=0; i < (2 * n_aug_ + 1); i++) {
      z_pred = z_pred + weights_(i) * Zsig_laser.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_laser, n_z_laser);
  S.fill(0.0);

  for (int i = 0; i < (2 * n_aug_ + 1); i++) {

    //residual
    VectorXd z_diff = Zsig_laser.col(i) - z_pred;

    S = S + (weights_(i) * z_diff * z_diff.transpose());

  }


  // Add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_laser, n_z_laser);
  R <<  (std_laspx_ * std_laspx_), 0,
        0, (std_laspy_ * std_laspy_);


  S = S + R;


  // Step - 2 : Use the actual measurement values to update the state mean and covariance matrix

  //create matrix for cross correlation Tc, between sigma points in state space and measurement space
  MatrixXd Tc = MatrixXd(n_x_, n_z_laser);
  Tc.fill(0.0);

  for (int i = 0; i < (2 * n_aug_ + 1); i++) {

    VectorXd z_diff = Zsig_laser.col(i) - z_pred;

    VectorXd x_diff = Xsig_pred.col(i) - x_;

    //angle normalization
    while (x_diff(3) > M_PI) {
      x_diff(3) -= (2. * M_PI);
    }

    while (x_diff(3) < -M_PI) {
      x_diff(3) += (2. * M_PI);
    }

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }


  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  VectorXd z_values = meas_package.raw_measurements_;

  // Diff between ground truth and mean predicted measurement values
  VectorXd z_diff = z_values - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + (K * z_diff);
  P_ = P_ - (K * S * K.transpose());

  nis_lidar = (z_diff.transpose() * S.inverse() * z_diff);
  std::cout<<"Lidar NIS: "<<nis_lidar<<endl;
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

  // Transform sigma points into measurement space (Radar)
  MatrixXd Zsig_radar = MatrixXd(n_z_radar, 2 * n_aug_ + 1);


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
  VectorXd z_pred = VectorXd(n_z_radar);
  z_pred.fill(0.0);

  for (int i=0; i < (2 * n_aug_ + 1); i++) {
      z_pred = z_pred + weights_(i) * Zsig_radar.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_radar, n_z_radar);
  S.fill(0.0);

  for (int i = 0; i < (2 * n_aug_ + 1); i++) {

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
  MatrixXd R = MatrixXd(n_z_radar, n_z_radar);
  R <<  (std_radr_ * std_radr_), 0, 0,
        0, (std_radphi_ * std_radphi_), 0,
        0, 0, (std_radrd_ * std_radrd_);


  S = S + R;


  // Step - 2 : Use the actual measurement values to update the state mean and covariance matrix

  //create matrix for cross correlation Tc, between sigma points in state space and measurement space
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    VectorXd z_diff = Zsig_radar.col(i) - z_pred;

    //angle normalization
    while (z_diff(1) > M_PI) {
      z_diff(1) -= (2. * M_PI);
    }

    while (z_diff(1) < -M_PI) {
      z_diff(1)+= (2. * M_PI);
    }

    VectorXd x_diff = Xsig_pred.col(i) - x_;

    //angle normalization
    while (x_diff(3) > M_PI) {
      x_diff(3) -= (2. * M_PI);
    }

    while (x_diff(3) < -M_PI) {
      x_diff(3) += (2. * M_PI);
    }

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }


  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  VectorXd z_values = meas_package.raw_measurements_;

  // Diff between ground truth and mean predicted measurement values
  VectorXd z_diff = z_values - z_pred;

  //angle normalization
  while (z_diff(1) > M_PI) {
    z_diff(1) -= (2. * M_PI);
  }

  while (z_diff(1) < -M_PI) {
    z_diff(1) += (2. * M_PI);
  }

  //update state mean and covariance matrix
  x_ = x_ + (K * z_diff);
  P_ = P_ - (K * S * K.transpose());

  nis_radar = (z_diff.transpose() * S.inverse() * z_diff);
  std::cout<<"Radar NIS: "<<nis_radar<<endl;
}
