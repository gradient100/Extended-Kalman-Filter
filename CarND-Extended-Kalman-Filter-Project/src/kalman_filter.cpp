#include "kalman_filter.h"
#include <cmath>
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
  /**
  TODO:
    * predict the state
  */
    x_ = F_ * x_;
    MatrixXd F_t = F_.transpose();
    P_ = F_ * P_ * F_t + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
    MatrixXd y = z - H_ * x_;
    MatrixXd H_t = H_.transpose();
    MatrixXd S = H_ * P_ * H_t + R_;
    MatrixXd S_i = S.inverse();
    MatrixXd K =  P_ * H_t * S_i;

    //new state
    x_ = x_ + (K * y);
    P_ = P_ - K * H_ * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
    // Normalize phi for z
    VectorXd z_copy = z;
    //float epsilon = 0;
    //float epsilon =0.005;
    float epsilon = 0.03;
    double phi = z[1];
    double pi = atan(1)*4;
    bool normalized = false;
    if (!(phi >= -1*pi-epsilon && phi <= pi+epsilon))
    {
      while (phi < -1*pi-epsilon)
          phi += 2*pi;
      while (phi > pi-epsilon)
          phi -= 2*pi; 
      normalized = true;
    }
    // if (phi >= pi)
    //   phi = pi;
    // if (phi <= -1*pi)
    //   phi = -1*pi;
  
    z_copy[1] = phi;

    // if (normalized)
    // {
    //   cout << " Normalized phi for z!" << endl;
    //   cout << "\t" << "Old phi : " << z[1] << endl;
    //   cout << "\t" << "New phi : " << z_copy[1] << endl;
    // }

    MatrixXd y = z_copy - h(x_);
    MatrixXd H_t = H_.transpose();
    MatrixXd S = H_ * P_ * H_t + R_;
    MatrixXd S_i = S.inverse();
    MatrixXd K =  P_ * H_t * S_i;

    //new state
    x_ = x_ + (K * y);
    P_ = P_- K * H_ * P_;
}

VectorXd KalmanFilter::h(const VectorXd &x) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  double px = x[0];
  double py = x[1];
  double vx = x[2];
  double vy = x[3];
  VectorXd z = VectorXd(3,1);
  z << 0.0, 0.0, 0.0;
  
  double rho = sqrt (px*px + py*py);
  if (fabs(rho) < 0.001) {
  //  cout << "rho = 0 when calculating h function";
    rho = 0.001;
    //return z; 
  }
  double phi;
  if (fabs(px) < 0.001)
    phi = atan2(py,0.001);
  else phi = atan2(py,px);
  // cout << "py = " << py << endl;
  // cout << "px = " << px << endl;
  // cout << "phi : " << phi << endl;
 
  double rho_dot = (px*vx + py*vy) / rho;
  // cout << "(px*vx + py*vy) = " << (px*vx + py*vy) << endl;
  // cout << "rho = " << rho << endl;
  // cout << "rho_dot = " << rho_dot << endl;

  
  z <<  rho, phi, rho_dot;
  return z;
}