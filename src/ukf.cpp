#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  is_initialized_ = false;
  time_us_ = 0;

  //Initialize state & augmented state dimension, lambda
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_.fill(0.0);

  //Weights

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if(!is_initialized_)
  {
    assert(meas_package.sensor_type_ == MeasurementPackage::LASER || meas_package.sensor_type_ == MeasurementPackage::RADAR);
    if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      double lidar_x = meas_package.raw_measurements_[0];
      double lidar_y = meas_package.raw_measurements_[1];
      x_ << lidar_x, lidar_y, 0, 0, 0;
      P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
            0, std_laspy_ * std_laspy_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      double ro = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double ro_dot = meas_package.raw_measurements_[2];
      x_ << ro * cos(phi), ro * sin(phi), 0, 0, 0;
      P_ << std_radr_ * std_radr_, 0, 0, 0, 0,
            0, std_radphi_ * std_radphi_, 0, 0, 0,
            0, 0, std_radrd_ * std_radrd_, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;   
    }
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }
  auto dt = (meas_package.timestamp_ - time_us_)/1000000.0; //(dt in seconds)

  //Prediction step
  Prediction(dt);

  //Update step using measurements
  if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    UpdateLidar(meas_package);
  else if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    UpdateRadar(meas_package);

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  //Create Augmented State Vector & covariance
  
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}