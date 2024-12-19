// Copyright (c) 2021 Franka Emika GmbH
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <riemannian_motion_policy/riemannian_motion_policy.hpp>
#include <cassert>
#include <cmath>
#include <exception>
#include <string>
#include <Eigen/Eigen>
#include <chrono>


using namespace std::chrono;

namespace {

template <class T, size_t N>
std::ostream& operator<<(std::ostream& ostream, const std::array<T, N>& array) {
  ostream << "[";
  std::copy(array.cbegin(), array.cend() - 1, std::ostream_iterator<T>(ostream, ","));
  std::copy(array.cend() - 1, array.cend(), std::ostream_iterator<T>(ostream));
  ostream << "]";
  return ostream;
}
}

namespace riemannian_motion_policy {

void RiemannianMotionPolicy::update_stiffness_and_references(){
  //target by filtering
  /** at the moment we do not use dynamic reconfigure and control the robot via D, K and T **/
  //K = filter_params_ * cartesian_stiffness_target_ + (1.0 - filter_params_) * K;
  //D = filter_params_ * cartesian_damping_target_ + (1.0 - filter_params_) * D;
  nullspace_stiffness_ = filter_params_ * nullspace_stiffness_target_ + (1.0 - filter_params_) * nullspace_stiffness_;
  //std::lock_guard<std::mutex> position_d_target_mutex_lock(position_and_orientation_d_target_mutex_);
  position_d_ = filter_params_ * position_d_target_ + (1.0 - filter_params_) * position_d_;
  
  orientation_d_ = orientation_d_.slerp(filter_params_, orientation_d_target_);
}

inline void pseudoInverse(const Eigen::MatrixXd& M_, Eigen::MatrixXd& M_pinv_, bool damped = true) {
  double lambda_ = damped ? 0.2 : 0.0;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(M_, Eigen::ComputeFullU | Eigen::ComputeFullV);   
  Eigen::JacobiSVD<Eigen::MatrixXd>::SingularValuesType sing_vals_ = svd.singularValues();
  Eigen::MatrixXd S_ = M_;  // copying the dimensions of M_, its content is not needed.
  S_.setZero();

  for (int i = 0; i < sing_vals_.size(); i++)
     S_(i, i) = (sing_vals_(i)) / (sing_vals_(i) * sing_vals_(i) + lambda_ * lambda_);

  M_pinv_ = Eigen::MatrixXd(svd.matrixV() * S_.transpose() * svd.matrixU().transpose());
}
//Calculate nearest point on sphere
Eigen::Vector3d RiemannianMotionPolicy::calculateNearestPointOnSphere(const Eigen::Vector3d& position,
                                                                      const Eigen::Vector3d& sphereCenter, 
                                                                      double radius) {
    // Calculate the vector from end-effector to sphere center
    Eigen::Vector3d vector = sphereCenter - position;

    // Calculate the magnitude of the vector
    double distanceToCenter = vector.norm();

    // Handle edge case where the end-effector is at the sphere center
    if (distanceToCenter < 1e-6) {
        std::cerr << "End-effector is at the sphere center. Returning sphere center as nearest point.\n";
        return sphereCenter;
    }

    // Normalize the vector to get the direction to the nearest point on the surface
    Eigen::Vector3d direction = vector / distanceToCenter;
    
    // Compute the nearest point on the sphere's surface to robot
    Eigen::Vector3d nearestPoint = vector - radius * direction;

    return nearestPoint;
}

//translate Jacobian to minimal distance point
Eigen::MatrixXd RiemannianMotionPolicy::calculate_J_obstacle(const Eigen::MatrixXd& link_jacobian_prev, const Eigen::MatrixXd& link_jacobian_next,const Eigen::VectorXd& points_link ,const Eigen::VectorXd& coord_link) {
  // Check for valid values
  Eigen::MatrixXd J_obstacle = Eigen::MatrixXd::Zero(6, 7);
  Eigen::MatrixXd J_p_obstacle = Eigen::MatrixXd::Zero(3, 7);
  Eigen::MatrixXd J_r_obstacle = Eigen::MatrixXd::Zero(3, 7);
  // Rotation part stays the same
  J_r_obstacle = link_jacobian_next.bottomRows(3);
  // Translation part
  J_p_obstacle = link_jacobian_prev.topRows(3);
  //make skewsymmetric matrix of coordlink
  Eigen::Vector3d r = points_link - coord_link;
  Eigen::Matrix3d r_skew;
  r_skew << 0, -r(2), r(1),
            r(2), 0, -r(0),
            -r(1), r(0), 0;
  J_p_obstacle = J_p_obstacle + r_skew * J_r_obstacle;
  J_obstacle.topRows(3) = J_p_obstacle;
  J_obstacle.bottomRows(3) = J_r_obstacle;
  return J_obstacle;
}


//RMP calculation for obstacle avoidance
Eigen::VectorXd RiemannianMotionPolicy::calculate_f_obstacle(const Eigen::VectorXd& d_obs, const Eigen::MatrixXd& Jp_obstacle) {
  
  Eigen::Matrix<double, 3, 1> nabla_d = Eigen::MatrixXd::Zero(3,1);
  double alpha_rep;
  double eta_rep = 3.8;
  double mu_rep = 0.5;
  double distance;
  Eigen::Matrix<double, 3, 1> f_repulsive = Eigen::MatrixXd::Zero(3,1);
  double alpha_damp;
  double eta_damp = 3.0;
  double mu_damp = 1.0;
  double epsilon = 0.001;
  Eigen::Matrix<double, 3, 1> P_obs = Eigen::MatrixXd::Zero(3,1);
  Eigen::Matrix<double, 3, 1> f_damping = Eigen::MatrixXd::Zero(3,1);
  Eigen::Matrix<double, 3, 1> f_obstacle = Eigen::MatrixXd::Zero(3,1);
  
  // Compute nabla_d
  nabla_d = d_obs/(d_obs.norm());
  distance = d_obs.norm();
  
  alpha_rep = eta_rep * std::exp(-distance / mu_rep);
  //alpha_rep = eta_rep * 1 / (std::pow(distance + mu_rep, 2));
  // Compute f_repulsive
  f_repulsive = alpha_rep * nabla_d.array();

  // Compute alpha_damp
  alpha_damp = eta_damp / ((distance / mu_damp) + epsilon);

  // Compute dot product for P_obs
  double dot_product = -(Jp_obstacle * dq_).transpose().dot(nabla_d);
  P_obs = std::max(0.0, dot_product) * nabla_d * nabla_d.transpose() * (Jp_obstacle * dq_);

  // Compute f_damping
  f_damping = -alpha_damp * P_obs.array();

  // Compute total f_obstacle
  f_obstacle = f_repulsive + f_damping;

  // Assign to f_obstacle_tilde
  Eigen::VectorXd f_obstacle_tilde = Eigen::VectorXd::Zero(6);
  f_obstacle_tilde.topRows(3) = f_obstacle;

  return f_obstacle_tilde;
}

Eigen::MatrixXd RiemannianMotionPolicy::calculate_A_obstacle(const Eigen::VectorXd& d_obs,
                                                             const Eigen::VectorXd& f_obstacle_tilde, double r_a) {
  double alpha_a = 1.0;
  double beta_x = 1.0;
  Eigen::Matrix3d H_obs = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d A_stretch = Eigen::Matrix3d::Identity();
  Eigen::Vector3d xsi = Eigen::Vector3d::Zero();
  Eigen::Vector3d f_obstacle = Eigen::Vector3d::Zero();
  Eigen::Matrix3d A_obs = Eigen::Matrix3d::Zero();
  Eigen::MatrixXd A_obs_tilde = Eigen::MatrixXd::Zero(6, 6);
  Eigen::Matrix3d identity_3 = Eigen::Matrix3d::Identity();

  f_obstacle = f_obstacle_tilde.topRows(3);

  // Check for valid values
  if (!d_obs.allFinite() || !f_obstacle.allFinite()) {
      throw std::runtime_error("d_obs or f_obstacle contains invalid values (NaN or inf).");
  }

  double c_1 = (-2.0 / r_a);
  double c_2 = 1.0 / std::pow(r_a, 2);
  double w_r;


  
  if (d_obs.norm() < r_a) {

    w_r = c_2 * d_obs.norm() * d_obs.norm() + c_1 * d_obs.norm() + 1.0;

  }
  else {
    w_r = 0.0;
  }

  double h_v = f_obstacle.norm() + log(1.0 + exp(-2 * alpha_a * f_obstacle.norm())) / alpha_a;
  xsi = f_obstacle / (h_v);

  A_stretch = xsi * xsi.transpose();
  H_obs = beta_x * A_stretch + (1.0 - beta_x) * identity_3;
  A_obs = w_r * H_obs ;
  
  A_obs_tilde.topLeftCorner(3, 3) = A_obs;

  return A_obs_tilde;
}


//RMP calculation of joint limit avoidance
void RiemannianMotionPolicy::rmp_joint_limit_avoidance(){
  //TODO: Implement the calculation of D_sigma fro joint limits
  //calculate sigma_u = 1/(1 + exp(-q))
  for (size_t i = 0; i < 7; ++i) {
    sigma_u(i) = 1/(1 + exp(-q_(i)));  
    //calculate alpha_u = 1/(1 + exp(-dq_ * c_alpha))
    alpha_u(i) = 1/(1 + exp(-dq_(i) * c_alpha));
    //calculate d_ii
    d_ii(i) = (q_upper_limit(i) - q_lower_limit(i)) * sigma_u(i) * (1 - sigma_u(i));
    //calculate d_ii_tilde
    d_ii_tilde(i) = sigma_u(i)*(alpha_u(i)*d_ii(i) + (1 - alpha_u(i))) + (1 - sigma_u(i))*((1- alpha_u(i)) * d_ii(i) + alpha_u(i));
    //calculate D_sigma
    D_sigma(i,i) = d_ii_tilde(i);
  }
  jacobian_tilde = jacobian * D_sigma;
  h_joint_limits = D_sigma.inverse() * (gamma_p * (q_0 - q_) - gamma_d * dq_);
}

void RiemannianMotionPolicy::calculateRMP_EE(const Eigen::MatrixXd& A_obs_tilde, 
                                                 const Eigen::VectorXd& f_obs_tilde) {
  
  Eigen::MatrixXd identity_3 = Eigen::MatrixXd::Identity(3, 3);
  Eigen::MatrixXd A_goal = Eigen::MatrixXd::Zero(6, 6);
  A_goal.topLeftCorner(3, 3) = identity_3;
  A_goal.bottomRightCorner(3, 3) = identity_3;
  A_tot = A_goal + A_obs_tilde;
    
  pseudoInverse(A_tot, A_tot_pinv);
    
  f_tot = A_tot_pinv*(A_goal * x_dd_des + A_obs_tilde * f_obs_tilde);
    
}



//get global joint acceleration for torque calculation
void RiemannianMotionPolicy::get_ddq(){

  Eigen::MatrixXd I_77 = Eigen::MatrixXd::Identity(7, 7);
  Eigen::MatrixXd I_66 = Eigen::MatrixXd::Identity(6, 6);
  Eigen::MatrixXd Link2_a = jacobian2tilde.transpose() * A_obs_tilde2 * jacobian2tilde;
  Eigen::MatrixXd Link3_a = jacobian3tilde.transpose() * A_obs_tilde3 * jacobian3tilde;
  Eigen::MatrixXd Link4_a = jacobian4tilde.transpose() * A_obs_tilde4 * jacobian4tilde;
  Eigen::MatrixXd Link5_a = jacobian5tilde.transpose() * A_obs_tilde5 * jacobian5tilde;
  Eigen::MatrixXd Link6_a = jacobian6tilde.transpose() * A_obs_tilde6 * jacobian6tilde;
  Eigen::MatrixXd Link7_a = jacobian7tilde.transpose() * A_obs_tilde7 * jacobian7tilde;
  Eigen::MatrixXd Hand_a = jacobianEEtilde.transpose() * A_obs_tildeEE* jacobianEEtilde;
  Eigen::MatrixXd Link2_b = jacobian2tilde.transpose() * A_obs_tilde2 * f_obs_tilde2;
  Eigen::MatrixXd Link3_b = jacobian3tilde.transpose() * A_obs_tilde3 * f_obs_tilde3;
  Eigen::MatrixXd Link4_b = jacobian4tilde.transpose() * A_obs_tilde4 * f_obs_tilde4;
  Eigen::MatrixXd Link5_b = jacobian5tilde.transpose() * A_obs_tilde5 * f_obs_tilde5;
  Eigen::MatrixXd Link6_b = jacobian6tilde.transpose() * A_obs_tilde6 * f_obs_tilde6;
  Eigen::MatrixXd Link7_b = jacobian7tilde.transpose() * A_obs_tilde7 * f_obs_tilde7;
  Eigen::MatrixXd Hand_b = jacobianEEtilde.transpose() * A_obs_tildeEE * f_obs_tildeEE;
  Eigen::MatrixXd A_total = jacobian_tilde.transpose()*I_66 *jacobian_tilde + Hand_a +Link2_a + Link3_a + Link4_a + Link5_a + Link6_a + Link7_a + lambda_RMP * I_77 * 1.5;
  Eigen::MatrixXd A_total_inv;
    
  A_total_inv = A_total.inverse();
  
  ddq_ = D_sigma * A_total_inv* (jacobian_tilde.transpose() * I_66 * x_dd_des + Hand_b +Link2_b + Link3_b + Link4_b + Link5_b + Link6_b + Link7_b + lambda_RMP * h_joint_limits * 1.5);
  
  //ddq_ = D_sigma * (jacobian_tilde.transpose()*A_tot *jacobian_tilde + lambda_RMP * I_77 ).inverse() * 
  //                  (jacobian_tilde.transpose() * A_tot * f_tot + lambda_RMP * h_joint_limits);
  
}

void RiemannianMotionPolicy::arrayToMatrix(const std::array<double,7>& inputArray, Eigen::Matrix<double,7,1>& resultMatrix)
{
 for(long unsigned int i = 0; i < 7; ++i){
     resultMatrix(i,0) = inputArray[i];
   }
}

void RiemannianMotionPolicy::arrayToMatrix(const std::array<double,6>& inputArray, Eigen::Matrix<double,6,1>& resultMatrix)
{
 for(long unsigned int i = 0; i < 6; ++i){
     resultMatrix(i,0) = inputArray[i];
   }
}

Eigen::Matrix<double, 7, 1> RiemannianMotionPolicy::saturateTorqueRate(
  const Eigen::Matrix<double, 7, 1>& tau_d_calculated,
  const Eigen::Matrix<double, 7, 1>& tau_J_d_M) {  
  Eigen::Matrix<double, 7, 1> tau_d_saturated{};
  for (size_t i = 0; i < 7; i++) {
  double difference = tau_d_calculated[i] - tau_J_d_M[i];
  tau_d_saturated[i] =
         tau_J_d_M[i] + std::max(std::min(difference, delta_tau_max_), -delta_tau_max_);
  }
  return tau_d_saturated;
}

controller_interface::InterfaceConfiguration
RiemannianMotionPolicy::command_interface_configuration() const {
  controller_interface::InterfaceConfiguration config;
  config.type = controller_interface::interface_configuration_type::INDIVIDUAL;
  for (int i = 1; i <= num_joints; ++i) {
    config.names.push_back(robot_name_ + "_joint" + std::to_string(i) + "/effort");
  }
  return config;
}


controller_interface::InterfaceConfiguration RiemannianMotionPolicy::state_interface_configuration()
  const {
  controller_interface::InterfaceConfiguration state_interfaces_config;
  state_interfaces_config.type = controller_interface::interface_configuration_type::INDIVIDUAL;

  for (int i = 1; i <= num_joints; ++i) {
    state_interfaces_config.names.push_back(robot_name_ + "_joint" + std::to_string(i) + "/position");
    state_interfaces_config.names.push_back(robot_name_ + "_joint" + std::to_string(i) + "/velocity");
  }

  for (const auto& franka_robot_model_name : franka_robot_model_->get_state_interface_names()) {
    state_interfaces_config.names.push_back(franka_robot_model_name);
    std::cout << franka_robot_model_name << std::endl;
  }

  const std::string full_interface_name = robot_name_ + "/" + state_interface_name_;

  return state_interfaces_config;
}


CallbackReturn RiemannianMotionPolicy::on_init() {
   UserInputServer input_server_obj(&position_d_target_, &rotation_d_target_, &K, &D, &T);
   std::thread input_thread(&UserInputServer::main, input_server_obj, 0, nullptr);
   input_thread.detach();
   return CallbackReturn::SUCCESS;
}


CallbackReturn RiemannianMotionPolicy::on_configure(const rclcpp_lifecycle::State& /*previous_state*/) {
  franka_robot_model_ = std::make_unique<franka_semantic_components::FrankaRobotModel>(
  franka_semantic_components::FrankaRobotModel(robot_name_ + "/" + k_robot_model_interface_name,
                                               robot_name_ + "/" + k_robot_state_interface_name));
                                               
  try {
    rclcpp::QoS qos_profile(1); // Depth of the message queue
    qos_profile.reliability(RMW_QOS_POLICY_RELIABILITY_RELIABLE);
    franka_state_subscriber = get_node()->create_subscription<franka_msgs::msg::FrankaRobotState>(
    "franka_robot_state_broadcaster/robot_state", qos_profile, 
    std::bind(&RiemannianMotionPolicy::topic_callback, this, std::placeholders::_1));
    std::cout << "Succesfully subscribed to robot_state_broadcaster" << std::endl;
  }

  catch (const std::exception& e) {
    fprintf(stderr,  "Exception thrown during publisher creation at configure stage with message : %s \n",e.what());
    return CallbackReturn::ERROR;
    }


  RCLCPP_DEBUG(get_node()->get_logger(), "configured successfully");
  return CallbackReturn::SUCCESS;
}


CallbackReturn RiemannianMotionPolicy::on_activate(
  const rclcpp_lifecycle::State& /*previous_state*/) {
  franka_robot_model_->assign_loaned_state_interfaces(state_interfaces_); 

  

  // Create the subscriber in the on_activate method
  desired_pose_sub = get_node()->create_subscription<geometry_msgs::msg::Pose>(
        "/riemannian_motion_policy/reference_pose", 
        10,  // Queue size
        std::bind(&RiemannianMotionPolicy::reference_pose_callback, this, std::placeholders::_1)
    );
  std::cout << "Succesfully subscribed to reference pose publisher" << std::endl;

  closest_point_sub_ = get_node()->create_subscription<messages_fr3::msg::ClosestPoint>(
        "/closest_point", 
        10,  // Queue size
        std::bind(&RiemannianMotionPolicy::closestPointCallback, this, std::placeholders::_1)
    );
  std::array<double, 16> initial_pose = franka_robot_model_->getPoseMatrix(franka::Frame::kEndEffector);
  Eigen::Affine3d initial_transform(Eigen::Matrix4d::Map(initial_pose.data()));
  position_d_ = initial_transform.translation();
  orientation_d_ = Eigen::Quaterniond(initial_transform.rotation());
  std::cout << "Completed Activation process" << std::endl;
  std::array<double, 7> gravity_force_vector_array = franka_robot_model_->getGravityForceVector();
  Eigen::Map<Eigen::Matrix<double, 7, 1>> gravity_force_vector(gravity_force_vector_array.data());
  tau_gravity = gravity_force_vector;
  return CallbackReturn::SUCCESS;
  
}


controller_interface::CallbackReturn RiemannianMotionPolicy::on_deactivate(
  const rclcpp_lifecycle::State& /*previous_state*/) {
  franka_robot_model_->release_interfaces();
  return CallbackReturn::SUCCESS;
}

std::array<double, 6> RiemannianMotionPolicy::convertToStdArray(const geometry_msgs::msg::WrenchStamped& wrench) {
    std::array<double, 6> result;
    result[0] = wrench.wrench.force.x;
    result[1] = wrench.wrench.force.y;
    result[2] = wrench.wrench.force.z;
    result[3] = wrench.wrench.torque.x;
    result[4] = wrench.wrench.torque.y;
    result[5] = wrench.wrench.torque.z;
    return result;
}

void RiemannianMotionPolicy::topic_callback(const std::shared_ptr<franka_msgs::msg::FrankaRobotState> msg) {
  // Existing handling of external forces
  O_F_ext_hat_K = convertToStdArray(msg->o_f_ext_hat_k);
  arrayToMatrix(O_F_ext_hat_K, O_F_ext_hat_K_M);
}


void RiemannianMotionPolicy::reference_pose_callback(const geometry_msgs::msg::Pose::SharedPtr msg)
{
    // Handle the incoming pose message
    std::cout << "received reference posistion as " <<  msg->position.x << ", " << msg->position.y << ", " << msg->position.z << std::endl;
    position_d_target_ << msg->position.x, msg->position.y,msg->position.z;
    orientation_d_target_.coeffs() << msg->orientation.x, msg->orientation.y, msg->orientation.z, msg->orientation.w;
    // You can add more processing logic here
}

void RiemannianMotionPolicy::closestPointCallback(const messages_fr3::msg::ClosestPoint::SharedPtr msg) {
    // Handle the incoming closest point message
    //std::cout << "received closest point as " <<  msg->x << ", " << msg->y << ", " << msg->z << std::endl;
    //std::cout << "closest frame is " << msg->frame << std::endl;
    d_obs2 << msg->frame2x, msg->frame2y, msg->frame2z;
    d_obs3 << msg->frame3x, msg->frame3y, msg->frame3z;
    d_obs4 << msg->frame4x, msg->frame4y, msg->frame4z;
    d_obs5 << msg->frame5x, msg->frame5y, msg->frame5z;
    d_obs6 << msg->frame6x, msg->frame6y, msg->frame6z;
    d_obs7 << msg->frame7x, msg->frame7y, msg->frame7z;
    d_obsEE << msg->frameeex, msg->frameeey, msg->frameeez;
    // Handle the coordinates of the closest point
    closest_point2 << msg->point2x, msg->point2y, msg->point2z;
    closest_point3 << msg->point3x, msg->point3y, msg->point3z;
    closest_point4 << msg->point4x, msg->point4y, msg->point4z;
    closest_point5 << msg->point5x, msg->point5y, msg->point5z;
    closest_point6 << msg->point6x, msg->point6y, msg->point6z;
    closest_point7 << msg->point7x, msg->point7y, msg->point7z;
    closest_pointEE << msg->pointeex, msg->pointeey, msg->pointeez;
    // Handle the Jacobian of the closest point
    jacobian_array2 = msg->jacobian2;
    jacobian_array3 = msg->jacobian3;
    jacobian_array4 = msg->jacobian4;
    jacobian_array5 = msg->jacobian5;
    jacobian_array6 = msg->jacobian6;
    jacobian_array7 = msg->jacobian7;
    jacobian_arrayEE = msg->jacobianee;
    //reshape to matrix (6x7) column major
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::ColMajor>> jacobian2obstacle(jacobian_array2.data());
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::ColMajor>> jacobian3obstacle(jacobian_array3.data());
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::ColMajor>> jacobian4obstacle(jacobian_array4.data());
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::ColMajor>> jacobian5obstacle(jacobian_array5.data());
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::ColMajor>> jacobian6obstacle(jacobian_array6.data());
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::ColMajor>> jacobian7obstacle(jacobian_array7.data());
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::ColMajor>> jacobianEEobstacle(jacobian_arrayEE.data());
    //assign to class variables
    jacobian2_obstacle = jacobian2obstacle;
    jacobian3_obstacle = jacobian3obstacle;
    jacobian4_obstacle = jacobian4obstacle;
    jacobian5_obstacle = jacobian5obstacle;
    jacobian6_obstacle = jacobian6obstacle;
    jacobian7_obstacle = jacobian7obstacle;
    jacobianEE_obstacle = jacobianEEobstacle;

}



void RiemannianMotionPolicy::jointStateCallback(const sensor_msgs::msg::JointState::SharedPtr msg) {
    // Check if the size of the effort vector is correct (should match the number of joints, e.g., 7 for Franka)
    if (msg->effort.size() == 7) {
        // Convert std::vector from the message into an Eigen matrix for tau_J
        for (size_t i = 0; i < 7; ++i) {
            tau_J(i) = msg->effort[i];  // Extract the measured joint torques
        }
    } else {
        RCLCPP_ERROR(get_node()->get_logger(), "JointState message has incorrect effort size");
    }
}


void RiemannianMotionPolicy::updateJointStates() {
  q_prev = q_;
  for (auto i = 0; i < num_joints; ++i) {
    const auto& position_interface = state_interfaces_.at(2 * i);
    const auto& velocity_interface = state_interfaces_.at(2 * i + 1);
    assert(position_interface.get_interface_name() == "position");
    assert(velocity_interface.get_interface_name() == "velocity");
    q_(i) = position_interface.get_value();
    dq_(i) = velocity_interface.get_value();
  }
}

controller_interface::return_type RiemannianMotionPolicy::update(const rclcpp::Time& /*time*/, const rclcpp::Duration& /*period*/) {  
  // if (outcounter == 0){
  // std::cout << "Enter 1 if you want to track a desired position or 2 if you want to use free floating with optionally shaped inertia" << std::endl;
  // std::cin >> mode_;
  // std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  // std::cout << "Mode selected" << std::endl;
  // while (mode_ != 1 && mode_ != 2){
  //   std::cout << "Invalid mode, try again" << std::endl;
  //   std::cin >> mode_;
  // }
  // }
  std::array<double, 49> mass = franka_robot_model_->getMassMatrix();
  std::array<double, 7> coriolis_array = franka_robot_model_->getCoriolisForceVector();
  std::array<double, 7> gravity_force_vector_array = franka_robot_model_->getGravityForceVector();
  //get pose matrix of each joint
  std::array<double, 16> pose2 = franka_robot_model_->getPoseMatrix(franka::Frame::kJoint2);
  std::array<double, 16> pose3 = franka_robot_model_->getPoseMatrix(franka::Frame::kJoint3);
  std::array<double, 16> pose4 = franka_robot_model_->getPoseMatrix(franka::Frame::kJoint4);
  std::array<double, 16> pose5 = franka_robot_model_->getPoseMatrix(franka::Frame::kJoint5);
  std::array<double, 16> pose6 = franka_robot_model_->getPoseMatrix(franka::Frame::kJoint6);
  std::array<double, 16> pose7 = franka_robot_model_->getPoseMatrix(franka::Frame::kJoint7);
  std::array<double, 16> poseEE = franka_robot_model_->getPoseMatrix(franka::Frame::kEndEffector);
  Eigen::Affine3d transform2(Eigen::Matrix4d::Map(pose2.data()));
  Eigen::Affine3d transform3(Eigen::Matrix4d::Map(pose3.data()));
  Eigen::Affine3d transform4(Eigen::Matrix4d::Map(pose4.data()));
  Eigen::Affine3d transform5(Eigen::Matrix4d::Map(pose5.data()));
  Eigen::Affine3d transform6(Eigen::Matrix4d::Map(pose6.data()));
  Eigen::Affine3d transform7(Eigen::Matrix4d::Map(pose7.data()));
  Eigen::Affine3d transformEE(Eigen::Matrix4d::Map(poseEE.data()));
  //get position of each joint
  joint2 = transform2.translation();
  joint3 = transform3.translation();
  joint4 = transform4.translation();
  joint5 = transform5.translation();
  joint6 = transform6.translation();
  joint7 = transform7.translation();
  jointEE = transformEE.translation();
  //jacobians
  

  Jp_obstacle2 = jacobian2_obstacle.topRows(3);
  Jp_obstacle3 = jacobian3_obstacle.topRows(3);
  Jp_obstacle4 = jacobian4_obstacle.topRows(3);
  Jp_obstacle5 = jacobian5_obstacle.topRows(3);
  Jp_obstacle6 = jacobian6_obstacle.topRows(3);
  Jp_obstacle7 = jacobian7_obstacle.topRows(3);
  Jp_obstacleEE = jacobianEE_obstacle.topRows(3);

  jacobian_array =  franka_robot_model_->getZeroJacobian(franka::Frame::kEndEffector);
  std::array<double, 16> pose = franka_robot_model_->getPoseMatrix(franka::Frame::kEndEffector);
  Eigen::Map<Eigen::Matrix<double, 7, 1>> coriolis(coriolis_array.data());
  Eigen::Map<Eigen::Matrix<double, 7, 1>> gravity_force_vector(gravity_force_vector_array.data());
  jacobian = Eigen::Map<Eigen::Matrix<double, 6, 7>> (jacobian_array.data());
  
  
  pseudoInverse(jacobian.transpose(), jacobian_transpose_pinv);
  pseudoInverse(jacobian, jacobian_pinv);
  Eigen::Map<Eigen::Matrix<double, 7, 7>> M(mass.data());
  Eigen::Affine3d transform(Eigen::Matrix4d::Map(pose.data()));
  Eigen::Vector3d position(transform.translation());
  Eigen::Quaterniond orientation(transform.rotation());
  orientation_d_target_ = Eigen::AngleAxisd(rotation_d_target_[0], Eigen::Vector3d::UnitX())
                        * Eigen::AngleAxisd(rotation_d_target_[1], Eigen::Vector3d::UnitY())
                        * Eigen::AngleAxisd(rotation_d_target_[2], Eigen::Vector3d::UnitZ());
  updateJointStates(); 

  
  error.head(3) << position - position_d_;

  if (orientation_d_.coeffs().dot(orientation.coeffs()) < 0.0) {
    orientation.coeffs() << -orientation.coeffs();
  }
  Eigen::Quaterniond error_quaternion(orientation.inverse() * orientation_d_);
  error.tail(3) << error_quaternion.x(), error_quaternion.y(), error_quaternion.z();
  error.tail(3) << -transform.rotation() * error.tail(3);

  error.head(3) << position - position_d_;

  //Calculate friction forces
  N = (Eigen::MatrixXd::Identity(7, 7) - jacobian_pinv * jacobian);
  Eigen::VectorXd  tau_nullspace(7), tau_d(7);
  pseudoInverse(jacobian.transpose(), jacobian_transpose_pinv);

  tau_nullspace <<  N * (nullspace_stiffness_ * config_control * (q_d_nullspace_ - q_) - //if config_control = true we control the whole robot configuration
                    (2.0 * sqrt(nullspace_stiffness_)) * dq_);  // if config control ) false we don't care about the joint position
  //d_obs1 = calculateNearestPointOnSphere(position, sphere_center, sphere_radius);
  //d_obs1 = d_obs_prev1 * 0.99 + d_obs1 * 0.01;
  Lambda = (jacobian * M.inverse() * jacobian.transpose()).inverse();
  x_dd_des = Lambda.inverse()*(-K_RMP * (error) - D_RMP * jacobian * dq_);
  f_obs_tildeEE = calculate_f_obstacle(d_obsEE, Jp_obstacleEE);
  A_obs_tildeEE = calculate_A_obstacle(d_obsEE, f_obs_tildeEE, 1.0);
  f_obs_tilde2 = calculate_f_obstacle(d_obs2, Jp_obstacle2);
  A_obs_tilde2 = calculate_A_obstacle(d_obs2, f_obs_tilde2, 0.1);
  f_obs_tilde3 = calculate_f_obstacle(d_obs3, Jp_obstacle3);
  A_obs_tilde3 = calculate_A_obstacle(d_obs3, f_obs_tilde3, 1.0);
  f_obs_tilde4 = calculate_f_obstacle(d_obs4, Jp_obstacle4);
  A_obs_tilde4 = calculate_A_obstacle(d_obs4, f_obs_tilde4, 1.0);
  f_obs_tilde5 = calculate_f_obstacle(d_obs5, Jp_obstacle5);
  A_obs_tilde5 = calculate_A_obstacle(d_obs5, f_obs_tilde5, 1.0);
  f_obs_tilde6 = calculate_f_obstacle(d_obs6, Jp_obstacle6);
  A_obs_tilde6 = calculate_A_obstacle(d_obs6, f_obs_tilde6, 1.0);
  f_obs_tilde7 = calculate_f_obstacle(d_obs7, Jp_obstacle7);
  A_obs_tilde7 = calculate_A_obstacle(d_obs7, f_obs_tilde7, 1.0);
  rmp_joint_limit_avoidance();
  jacobian2tilde = jacobian2_obstacle * D_sigma;
  jacobian3tilde = jacobian3_obstacle * D_sigma;
  jacobian4tilde = jacobian4_obstacle * D_sigma;
  jacobian5tilde = jacobian5_obstacle * D_sigma;
  jacobian6tilde = jacobian6_obstacle * D_sigma;
  jacobian7tilde = jacobian7_obstacle * D_sigma;
  jacobianEEtilde = jacobianEE_obstacle * D_sigma;
  get_ddq();
  
  // Calculate the desired torque
  tau_RMP = M * ddq_;
  // Calculate friction torques
  calculate_tau_friction();
  calculate_tau_gravity(coriolis, gravity_force_vector, jacobian);
  //tau_gravity_error = tau_gravity - gravity_force_vector;

  auto tau_d_placeholder = tau_RMP + coriolis + tau_friction; //add nullspace, friction, gravity and coriolis components to desired torque
  tau_d << tau_d_placeholder;
  tau_d << saturateTorqueRate(tau_d, tau_J_d_M);  // Saturate torque rate to avoid discontinuities
  tau_J_d_M = tau_d;

  // Step 5: Implement a logger that logs every 5 seconds
  static steady_clock::time_point last_log_time = steady_clock::now();
  steady_clock::time_point current_time = steady_clock::now();

  // Check if 5 seconds have passed since the last log
  if (duration_cast<seconds>(current_time - last_log_time).count() >= 0.1) {
    last_log_time = current_time;  // Reset the last log time
        
  }

  for (size_t i = 0; i < 7; ++i) {
    command_interfaces_[i].set_value(tau_d(i));
  }
  
  if (outcounter % 1000/update_frequency == 0){
    std::cout << "x_dd_des" << std::endl;
    std::cout << x_dd_des << std::endl;
    std::cout << "error_pose" << std::endl;
    std::cout << error << std::endl;
    std::cout << "beta_orientation" << std::endl;
    std::cout << beta_orientation << std::endl;
    std::cout << "dq_" << std::endl;
    std::cout << dq_ << std::endl;
    std::cout << "ddq_" << std::endl;
    std::cout << ddq_ << std::endl;
    std::cout << "tau_RMP" << std::endl;
    std::cout << tau_RMP << std::endl;
    std::cout << "d_obs2" << std::endl;
    std::cout << d_obs2 << std::endl;
    std::cout << "d_obs3" << std::endl;
    std::cout << d_obs3 << std::endl;
    std::cout << "d_obs4" << std::endl;
    std::cout << d_obs4 << std::endl;
    std::cout << "d_obs5" << std::endl;
    std::cout << d_obs5 << std::endl;
    std::cout << "d_obs6" << std::endl;
    std::cout << d_obs6 << std::endl;
    std::cout << "d_obs7" << std::endl;
    std::cout << d_obs7 << std::endl;
    std::cout << "d_obsEE" << std::endl;
    std::cout << d_obsEE << std::endl;
    std::cout << "jointee" << std::endl;
    std::cout << jointEE << std::endl;
    std::cout << "jacobian_obstacle" << std::endl;
    std::cout << jacobianEE_obstacle << std::endl;

  }
  outcounter++;
  update_stiffness_and_references();
  return controller_interface::return_type::OK;
}
}

// namespace riemannian_motion_policy
#include "pluginlib/class_list_macros.hpp"
// NOLINTNEXTLINE
PLUGINLIB_EXPORT_CLASS(riemannian_motion_policy::RiemannianMotionPolicy,
                       controller_interface::ControllerInterface)