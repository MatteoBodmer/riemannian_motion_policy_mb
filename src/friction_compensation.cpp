#include <riemannian_motion_policy/riemannian_motion_policy.hpp>

namespace riemannian_motion_policy{
    //Calculates the friction forces acting on the robot's joints depending on joint rotational speed. 
    //Exerts torque up to a certain empirically detected static friction threshold. 
    //TODO: Afterwards, goes into the viscous domain and follows a linear raise depending on empiric parameters
    void RiemannianMotionPolicy::calculate_tau_friction(){
        if (friction_){
            double alpha = 0.01;//constant for exponential filter in relation to static friction moment        
            dq_filtered = alpha* dq_ + (1 - alpha) * dq_filtered; //Filtering dq of every joint
            tau_impedance_filtered = alpha* tau_impedance + (1 - alpha) * tau_impedance_filtered; //Filtering tau_admittance
            //Creating and filtering a "fake" tau_admittance with own weights, optimized for friction compensation (else friction compensation would get stronger with higher stiffnesses)
            tau_threshold = jacobian.transpose() * Sm * (-alpha * (D_friction*(jacobian*dq_) + K_friction * error)) + (1 - alpha) * tau_threshold;
            //Creating "fake" dq, that acts only in the impedance-space, else dq in the nullspace also gets compensated, which we do not want due to null-space movement
            Eigen::Matrix<double, 7, 1> dq_compensated = dq_filtered - N * dq_filtered;

            //Calculation of friction force according to Bachelor Thesis: https://polybox.ethz.ch/index.php/s/iYj8ALPijKTAC2z?path=%2FFriction%20compensation
            f = beta.cwiseProduct(dq_compensated) + offset_friction;
            g = coulomb_friction;
            g(4) = (coulomb_friction(4) + (static_friction(4) - coulomb_friction(4)) * exp(-1 * std::abs(dq_compensated(4)/dq_s(4))));
            g(6) = (coulomb_friction(6) + (static_friction(6) - coulomb_friction(6)) * exp(-1 * std::abs(dq_compensated(6)/dq_s(6))));
           
            dz = dq_compensated.array() - dq_compensated.array().abs() / g.array() * sigma_0.array() * z.array() + 0.025* tau_threshold.array()/*(jacobian.transpose() * K * error).array()*/;
            /** 
            std::cout << "g is " << g.transpose() << std::endl;
            std::cout << "dz is " << dz.transpose() << std::endl;
            std::cout << "z is " << dz.transpose() << std::endl;
            std::cout << "dq_compensated is " << dq_compensated.transpose() << std::endl;
            std::cout << "tau thresh is " << tau_threshold.transpose() << std::endl;
            */
            dz(6) -= 0.02*tau_threshold(6); // dz is nana
            z = 0.001 * dz + z;
            tau_friction = sigma_0.array() * z.array() + 100 * sigma_1.array() * dz.array() + f.array();  
            
        }

        else
        {
            tau_friction.setZero();

        }
        
    }
}