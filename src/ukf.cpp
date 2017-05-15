#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
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
  std_a_ = 1.3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.515;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
   
    n_x_ = 5;
    n_aug_ = 7;
    Xsig_pred_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    lambda_ = 3 - n_aug_;
    weights_ = VectorXd(2*n_aug_+1);
    
    P_ <<   1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
    
    //initialize weights
    weights_(0) = lambda_/(lambda_+n_aug_);
    for(int i=1;i<2*n_aug_+1;i++)
    {
        double weight = 0.5/(n_aug_+lambda_);
        weights_(i) = weight;
    }
    NIS_laser_ = 0;
    NIS_radar_ = 0;
    
    previous_timestamp_ = 0;
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
        x_ = VectorXd(5);
        
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            float rho = meas_package.raw_measurements_(0);
            float phi = meas_package.raw_measurements_(1);
            
            float px = rho * cos(phi);
            float py = rho * sin(phi);
            
            x_ << px, py, 0, 0, 0;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            float px = meas_package.raw_measurements_(0);
            float py = meas_package.raw_measurements_(1);
            x_ << px, py, 0, 0, 0;
        }
        
      previous_timestamp_ = meas_package.timestamp_;
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }
    
    double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;

    
    Prediction(dt);
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
    {
        UpdateRadar(meas_package);
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
    {
        UpdateLidar(meas_package);
    }
    previous_timestamp_ = meas_package.timestamp_;    
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    
    //Generate Sigma Points
    
    //create augmented mean vector
    VectorXd x_aug = VectorXd(7);
    
    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(7,7);
    
    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
    
    //create augmented mean state
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;
    
    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5, 5) = P_;
    P_aug(5,5) = std_a_ * std_a_;
    P_aug(6,6) = std_yawdd_ * std_yawdd_;
    
    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();
    
    //create augmented sigma points
    Xsig_aug.col(0) = x_aug;
    for(int i = 0; i < n_aug_;i++)
    {
        Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }
    
    
    
    //predict sigma points
    // predict sigma points
    for (int i=0; i<2*n_aug_+1; i++) {
        
        double px = Xsig_aug(0,i);
        double py = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawRate = Xsig_aug(4,i);
        double v_acc_noise = Xsig_aug(5,i);
        double yaw_acc_noise = Xsig_aug(6,i);
        
        double px_pred;
        double py_pred;
        double v_pred;
        double yaw_pred;
        double yawRate_pred;
        
        double delta_t_sq = delta_t*delta_t;
        
        yaw_pred = yaw + yawRate*delta_t;
        
        if (fabs(yawRate) < 0.001) {
            px_pred = px + (v*cos(yaw)*delta_t);
            py_pred = py + (v*sin(yaw)*delta_t);
        } else {
            px_pred = px + ((v/yawRate) * (sin(yaw_pred)-sin(yaw)));
            py_pred = py + ((v/yawRate) * (-cos(yaw_pred)+cos(yaw)));
        }
        
        // noise
        px_pred  += 0.5 * delta_t_sq * cos(yaw) * v_acc_noise;
        py_pred  += 0.5 * delta_t_sq * sin(yaw) * v_acc_noise;
        yaw_pred += 0.5 * delta_t_sq * yaw_acc_noise;
        yawRate_pred = yawRate + delta_t * yaw_acc_noise;
        v_pred = v + delta_t * v_acc_noise;
        
        // assign predictions
        Xsig_pred_(0,i) = px_pred;
        Xsig_pred_(1,i) = py_pred;
        Xsig_pred_(2,i) = v_pred;
        Xsig_pred_(3,i) = yaw_pred;
        Xsig_pred_(4,i) = yawRate_pred;
    }
    
    // predict state
    x_.fill(0.0);
    for (int i=0; i<2*n_aug_+1; i++) {
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }
    
    //predict state covraince matrix
    //create covariance matrix for prediction
    P_.fill(0.0);
    for(int i=0;i<2*n_aug_+1;i++)
    {
        VectorXd x_diff = Xsig_pred_.col(i)  - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
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
    
    int n_z = 2;
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    
    for (int i=0; i<2*n_aug_+1;i++)
    {
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        
        // measurement model
        Zsig(0,i) = p_x;                        //r
        Zsig(1,i) = p_y;                                 //phi
    }
    
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R << std_laspx_*std_laspx_, 0,  0, std_laspy_ * std_laspy_;
    S = S + R;
    
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    VectorXd z = VectorXd(n_z);
    z << meas_package.raw_measurements_;
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    
    //residual
    VectorXd z_diff = z - z_pred;
    
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    
    //NIS radar calculation
    
    NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
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
   //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    
    for (int i=0; i<2*n_aug_+1;i++)
    {
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
        
        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
        Zsig(1,i) = atan2(p_y,p_x);                                 //phi
        Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
    
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    
    
    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0,std_radrd_*std_radrd_;
    S = S + R;
    
    
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    VectorXd z = VectorXd(n_z);
    z << meas_package.raw_measurements_;
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    
    //residual
    VectorXd z_diff = z - z_pred;
    
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    
    //NIS radar calculation
    
    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
