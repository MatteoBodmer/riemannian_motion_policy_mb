#include <riemannian_motion_policy/riemannian_motion_policy.hpp>

namespace riemannian_motion_policy{
    //Calculates the gravity compensation torque.

    void RiemannianMotionPolicy::calculate_residual_torque(const Eigen::Map<Eigen::Matrix<double, 7, 1>>& coriolis, const Eigen::Map<Eigen::Matrix<double, 7, 1>>& gravity_force_vector) {
        // Use previous timestep to calculate M_dot
        K_0.diagonal() << 80,80,80,80,80,40,40;
        
        // Low-pass filter for dq_

        Eigen::Matrix<double, 7, 7> M_dot = (M - M_prev) / dt;

    
        I_tau += (tau_J - gravity_force_vector + M_dot * dq_  - coriolis + residual_prev) * dt; 

                   
        // Compute the residual
        residual = K_0 * (M * dq_ - I_tau);
        // Update previous values
        residual_prev = residual;
        M_prev = M;
    }


    // Function to filter the residual torque and extract gravity torque using RTOB
    void RiemannianMotionPolicy::calculateRTOB() {
        // Apply low-pass filter to the residual torque (recursive torque observer)
        tau_gravity_filtered = filter_gain * residual + (1.0 - filter_gain) * tau_gravity_filtered_prev;

        // Update the previous filtered value for the next iteration
        tau_gravity_filtered_prev = tau_gravity_filtered;
    }


    // Function to calculate gravity torques for all joints
    /*void RiemannianMotionPolicy::calculate_gravity_torques() {
        // Initialize the intercepts for each joint
        Eigen::Matrix<double, 7, 1> intercepts;
        intercepts << -0.0004, -32.9103, -0.3486, 13.9538, 0.8066, 2.3977, -0.0120;

    // Initialize the coefficients for each joint (7 joints, 7 coefficients per joint for q1 to q7)
        Eigen::Matrix<double, 7, 7> coefficients;
        coefficients << 
            -0.0004,  0.0006, -0.0004, -0.0008,  0.0002, -0.0008,  0.0005,  // Joint 1
            -10.4280, -76.7979,  2.7157, 41.0974, -14.1278, 36.7995, 10.1022,  // Joint 2
            -0.8839,  4.7166, -3.1532, -1.0539, -20.3315, -0.6610,  0.5505,  // Joint 3
            1.7845,  4.7020,  4.6835, -6.1887, -9.5141, -0.3636, -1.8354,  // Joint 4
            0.0992,  0.1803,  0.1069, -0.1230,  0.0692, -0.1229, -0.0146,  // Joint 5
            0.1628, -0.8906,  0.1808,  0.7837, -0.4687,  0.7331, -0.1956,  // Joint 6
            -0.0114,  0.0225, -0.0131, -0.0154,  0.0204, -0.0127,  0.0113;  // Joint 7



        // Calculate the gravity torques: intercepts + coefficients * q_
        tau_gravity = intercepts + coefficients * q_;
    }*/

    void RiemannianMotionPolicy::calculate_tau_gravity(const Eigen::Map<Eigen::Matrix<double, 7, 1>>& coriolis, const Eigen::Map<Eigen::Matrix<double, 7, 1>>& gravity_force_vector, const Eigen::Matrix<double, 6, 7>& jacobian) {
        if (gravity_){
                   
                // Step 1: Calculate disturbance torques using the 1st Momentum Observer Residual
                calculate_residual_torque(coriolis, gravity_force_vector);  // Calculate the residual torque
                // Step 2: Apply RTOB-specific filtering to extract the gravity component from the disturbance    
                calculateRTOB();  // Filtering in RTOB
                // Step 3: Calculate gravity torques using the model
                tau_gravity = tau_gravity_filtered;

               // Final gravity torque 
            
               
            
        }

        else
        {
            tau_gravity_error.setZero();
        }
        
    }




    

}

