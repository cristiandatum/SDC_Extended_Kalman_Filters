#include "kalman_filter.h"
#include <iostream>
#define PI 3.14159265

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {

    // Time Update: Extrapolate the state (new estimate)
    x_ = (F_ * x_); // + u; removed for external motion

    // Time update: Extrapolate uncertainty (new uncertainty)
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

    // Measurement Update: Kalman Gain equation:
    MatrixXd K = P_ * H_.transpose() * (H_ * P_ * H_.transpose() + R_).inverse();

    // Measurement Update: update estimate with measurement:
    x_ = x_ + K * (z - (H_ * x_)); 

    // Determine Identity matrix required size
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);

    // Measurement Update: update the estimate with uncertainty
    P_ = (I - K * H_) * P_ * (I - (K * H_)).transpose() + K * R_ * (K.transpose());

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

    Tools tools;
    MatrixXd Hj_ = tools.CalculateJacobian(x_);

    float px = x_[0];
    float py = x_[1];
    float vx = x_[2];
    float vy = x_[3];
    float rho = sqrt(px*px + py*py);
    float phi = atan2f(py, px);
    float rho_dot = (px*vx + py*vy)/rho;

    VectorXd hx_(3);
    VectorXd y_(3);

    hx_ << rho, phi, rho_dot;

    // Measurement Update: Extended Kalman Gain equation:
    MatrixXd K = P_ * Hj_.transpose() * (Hj_ * P_ * Hj_.transpose() + R_).inverse();

    // Measurement Update: update estimate with measurement:
    y_ = z - hx_;
    
    // Check if difference between angles is greater than PI, or lesser than -PI
    if (y_[1]>PI){
        y_[1]=y_[1]-PI;
    } 
    if (y_[1]<-PI){
        y_[1]=y_[1]+PI;
    } 
    x_ = x_ + K * y_;


    // Determine Identity matrix required size
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);

    // Measurement Update: update the estimate with uncertainty
    P_ = (I - K * Hj_) * P_ * (I - (K * Hj_)).transpose() + K * R_ * (K.transpose());

}
