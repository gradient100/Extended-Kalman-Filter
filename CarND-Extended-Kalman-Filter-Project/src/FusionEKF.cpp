#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  timeStep = 1;


  // initializing measurement covariance matrices
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0, 0.0225;
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ = MatrixXd(2, 4);
  H_laser_ <<  1, 0, 0, 0,
               0, 1, 0, 0;
  Hj_ = MatrixXd(3, 4);

  ekf_.x_ = VectorXd(4);
  ekf_.F_ = MatrixXd(4,4);  
  ekf_.P_ = MatrixXd(4,4);
  ekf_.Q_ = MatrixXd(4, 4);


  // ekf_.H_ and ekf_.R_ initialization in measurement function (size known when sensor known)
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  //if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) 
  //if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  
  {
    if (!is_initialized_) {
      /**
      TODO:
        * Initialize the state ekf_.x_ with the first measurement.
        * Create the covariance matrix.
        * Remember: you'll need to convert radar from polar to cartesian coordinates.
      */
      // first measurement
      cout << "EKF: " << endl;
      //ekf_.x_ = VectorXd(4);
      //ekf_.x_ << 1, 1, 1, 1;


      ekf_.P_ <<  1, 0, 0, 0,
                  0, 1, 0, 0, 
                  0, 0, 1000, 0, 
                  0, 0, 0, 1000;

      cout << " Time Step = " << timeStep++ << endl;
      if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        /**
        Convert radar from polar to cartesian coordinates and initialize state.
        */
        
        cout << "***** RADAR (initializing) ******" << endl;
        double rho = measurement_pack.raw_measurements_[0];
        double phi = measurement_pack.raw_measurements_[1];
        double pi = atan(1)*4;

        while (!(phi >= -1*pi && phi <= pi))
        {
          while (phi < -1*pi)
              phi += 2*pi;
          while (phi > pi)
              phi -= 2*pi; 
        }
        ekf_.x_ << rho * cos(phi), rho * sin(phi), 0, 0;

        cout << "phi : " << phi << endl;
        cout << "x_ initialized = " << endl << ekf_.x_ << endl;
        cout << "P_ initialized = " << endl << ekf_.P_ << endl;

      }

      else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        /**
        Initialize state.
        */
        
        cout << "***** LASER (initializing) *****" << endl;
        ekf_.x_ << measurement_pack.raw_measurements_[0], 
                    measurement_pack.raw_measurements_[1], 0, 0;
        cout << "x_ initialized = " << endl << ekf_.x_ << endl;
        cout << "P_ initialized = " << endl << ekf_.P_ << endl;
      
      }

      // done initializing, no need to predict or update
      is_initialized_ = true;
      return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    /**
     TODO:
       * Update the state transition matrix F according to the new elapsed time.
        - Time is measured in seconds.
       * Update the process noise covariance matrix.
       * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */
    //Modify the F matrix so that the time is integrated
    
    if (previous_timestamp_==0)
        previous_timestamp_ = measurement_pack.timestamp_;
    double dt = (double)(measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
    double dt_2 = dt * dt;
    double dt_3 = dt_2 * dt;
    double dt_4 = dt_3 * dt;
    
    previous_timestamp_ = measurement_pack.timestamp_;

    ekf_.F_ << 1, 0, dt, 0,
              0, 1, 0, dt,
              0, 0, 1, 0,
              0, 0, 0, 1;


    double noise_ax = 9;
    double noise_ay = 9;
    //set the process covariance matrix Q
    ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
           0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
           dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
           0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

    cout << endl << "**********************************************************************************" <<endl;
    cout << " Time Step = " << timeStep++ << endl;
    cout << endl << "dt : " << dt << endl;

    ekf_.Predict();

    cout << "x_ predicted = " << endl << ekf_.x_ << endl;
    cout << "P_ predicted = " << endl << ekf_.P_ << endl;

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
     TODO:
       * Use the sensor type to perform the update step.
       * Update the state and covariance matrices.
     */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Radar updates
      cout << endl << "***** RADAR *****" << endl;

      ekf_.R_ = R_radar_;
      Hj_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.H_ = Hj_;

      

      ekf_.UpdateEKF(measurement_pack.raw_measurements_);

    } 
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Laser updates
      cout << endl << "***** LASER *****" << endl;

      ekf_.R_ = R_laser_;
      ekf_.H_ = H_laser_;


      ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ updated = " << endl << ekf_.x_ << endl;
    cout << "P_ updated = " << endl << ekf_.P_ << endl;
  }
}
