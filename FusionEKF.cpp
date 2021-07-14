#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // Measurement covariance matrix - laser
    R_laser_ = MatrixXd(2, 2);
    R_laser_ << 0.0225, 0,
                0, 0.0225;

    // Measurement covariance matrix - radar
    R_radar_ = MatrixXd(3, 3);
    R_radar_ << 0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;

    // H_laser measurement matrix:
    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
                0, 1,  0, 0;
  
    // H_j radar measurement matrix:
    Hj_ = MatrixXd(3, 4);
    Hj_ << 1 /sqrt(2), 1/sqrt(2), 0, 0,
           -0.5, -0.5, 0, 0,
           0, 0, 1, 1;

    // State covariance matrix P
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;

    // Process covariance matrix Q
    ekf_.Q_ = MatrixXd(4,4);
    ekf_.Q_ << 0, 0, 0, 0,
               0, 0, 0, 0,
               0, 0, 0, 0,
               0, 0, 0, 0; 

    // Acceleration noise components
    float noise_ax = 9;
    float noise_ay = 9;

    // Initial transition matrix
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;

    // Initial state measurement (px, py, vx, vy)
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 0, 0;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    /**
     * Initialization
     */
    if (!is_initialized_) {

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // TODO: Convert radar from polar to cartesian coordinates 
        //         and initialize state.
            ekf_.x_ << measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[0]),
                       measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]),
                       0, 
                       0;
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        // TODO: Initialize state.
            ekf_.x_ << measurement_pack.raw_measurements_[0], 
                       measurement_pack.raw_measurements_[1], 
                       0, 
                       0;
        }
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
    /*** Time Update ***/
    // compute the time elapsed between the current and previous measurements
    // dt - expressed in seconds
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    ekf_.F_ << 1, 0, dt, 0,
                0, 1, 0, dt,
                0, 0, 1, 0,
                0, 0, 0, 1;

    ekf_.Q_ << 0.25 * dt * dt * dt * dt * noise_ax, 0, 0.5 * dt * dt * dt * noise_ax, 0,
                0, 0.25 * dt * dt * dt * dt * noise_ay, 0 , 0.5 * dt * dt * dt * noise_ay,
                0.5 * dt * dt * dt * noise_ax, 0, dt * dt * noise_ax, 0,
                0, 0.5 * dt * dt * dt * noise_ay, 0, dt * dt * noise_ay; 

    ekf_.Predict();

    /*** Measurement Update ***/

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
        ekf_.R_ = R_radar_;
        ekf_.H_ = Hj_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
        cout << "updating radar"<< endl;
  } else {
    // Laser updates
        ekf_.R_ = R_laser_;
        ekf_.H_ = H_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
        cout << "updating laser"<< endl;

  }

  // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
