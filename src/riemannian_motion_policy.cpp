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
#include <yaml-cpp/yaml.h>>


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

//RMP calculation for obstacle avoidance
Eigen::VectorXd RiemannianMotionPolicy::calculate_f_obstacle(const Eigen::VectorXd& d_obs, const Eigen::MatrixXd& Jp_obstacle) {
  
  Eigen::Matrix<double, 3, 1> nabla_d = Eigen::MatrixXd::Zero(3,1);
  double alpha_rep; //repulsive potential scalar
  double alpha_damp; //repulsive damping scalar
  double distance; //from given point to nearest obstacle point
  Eigen::Matrix<double, 3, 1> f_repulsive = Eigen::MatrixXd::Zero(3,1);
  Eigen::Matrix<double, 3, 1> P_obs = Eigen::MatrixXd::Zero(3,1);
  Eigen::Matrix<double, 3, 1> f_damping = Eigen::MatrixXd::Zero(3,1);
  Eigen::Matrix<double, 3, 1> f_obstacle = Eigen::MatrixXd::Zero(3,1);
  Eigen::Matrix<double, 3, 1> w = Eigen::MatrixXd::Zero(3,1); // task space linear velocity
  
  // Compute nabla_d
  nabla_d = d_obs/(std::max(d_obs.norm(),0.001));
  distance = d_obs.norm();
  w = Jp_obstacle * dq_;
  
  // Compute f_repulsive
  alpha_rep = eta_rep * std::exp(-distance / mu_rep);
  f_repulsive = alpha_rep * nabla_d.array();

  // Compute alpha_damp
  alpha_damp = eta_damp / ((distance / mu_damp) + epsilon);
  // stretching matrix (only consider velocities when moving towards the obstacle)
  P_obs = std::max(0.0, -w.transpose().dot(nabla_d)) * nabla_d * nabla_d.transpose() * w;
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
                                                             const Eigen::VectorXd& f_obstacle_tilde, double r_a,
                                                             const Eigen::MatrixXd& Jp_obstacle) {
  double alpha_a = 1.0; // tuning paramter for directional stretching
  double beta_x = 1.0; // tuning parameter wight for directional stretching vs Identity
  
  Eigen::Matrix3d H_obs = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d A_stretch = Eigen::Matrix3d::Identity();
  Eigen::Vector3d xi = Eigen::Vector3d::Zero();
  Eigen::Vector3d f_obstacle = Eigen::Vector3d::Zero();
  Eigen::Matrix3d A_obs = Eigen::Matrix3d::Zero();
  Eigen::MatrixXd A_obs_tilde = Eigen::MatrixXd::Zero(6, 6);
  Eigen::Matrix3d identity_3 = Eigen::Matrix3d::Identity();
  Eigen::Matrix<double, 3, 1> nabla_d = Eigen::MatrixXd::Zero(3,1);
  Eigen::Matrix<double, 3, 1> w = Eigen::MatrixXd::Zero(3,1); // task space linear velocity
  // define interpolation polynomial (omega)
  double c_1 = (-2.0 / r_a);
  double c_2 = 1.0 / std::pow(r_a, 2);
  double w_r;

  // define linear components of repulsion force
  f_obstacle = f_obstacle_tilde.topRows(3);
  // define repulsion direction and linear velocity
  nabla_d = d_obs/(std::max(d_obs.norm(), 0.001));
  w = Jp_obstacle * dq_;

  // Check for valid values
  if (!d_obs.allFinite() || !f_obstacle.allFinite()) {
      throw std::runtime_error("d_obs or f_obstacle contains invalid values (NaN or inf).");
  }
  if (d_obs.norm() < r_a) {
    w_r = c_2 * d_obs.norm() * d_obs.norm() + c_1 * d_obs.norm() + 1.0;

  }
  else {
    w_r = 0.0;
  }

  double h_v = f_obstacle.norm() + 1/alpha_a * log(1.0 + exp(-2 * alpha_a * f_obstacle.norm()));
  xi = f_obstacle / (h_v);
  A_stretch = xi * xi.transpose();

  H_obs = beta_x * A_stretch + (1.0 - beta_x) * identity_3;
  A_obs = w_r * H_obs ;
  A_obs_tilde.topLeftCorner(3, 3) = A_obs;

  return weight_obstacle * A_obs_tilde;
}

//RMP for target attraction/axis
Eigen::MatrixXd RiemannianMotionPolicy::calculate_target_attraction(const Eigen::VectorXd& error, const Eigen::MatrixXd& jacobian) {
  //declarations
  Eigen::VectorXd f_attract = Eigen::VectorXd::Zero(6);  
  Eigen::MatrixXd A_attract = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd j_translational = jacobian.topRows(3);
  Eigen::MatrixXd j_rotational = jacobian.bottomRows(3);
  Eigen::Vector3d error_position = error.head(3);
  Eigen::Vector3d error_orientation = error.tail(3);
  Eigen::Matrix3d M_far = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d M_near = Eigen::Matrix3d::Zero();
  //target metric
  Eigen::Matrix3d A_position = Eigen::Matrix3d::Zero();
  double alpha = (1 - alpha_min) *exp((-1 * error_position.squaredNorm()) / (2*sigma_a)) + alpha_min;
  double beta = exp((-1 * error_position.squaredNorm()) / (2*sigma_b));
  M_near = Eigen::Matrix3d::Identity();
  M_far = 0/(error_position.squaredNorm()) * error_position * error_position.transpose();
  A_position = (beta * b + (1 - beta)) * (alpha * M_near + (1 - alpha) * M_far);
  A_attract.topLeftCorner(3, 3) = A_position;
  
  //axis metric
  Eigen::Matrix3d A_orientation = Eigen::Matrix3d::Zero();
  double beta_axis = exp((-1 * std::pow(error_position.norm(), 2)) / (2*sigma_o));
  A_orientation = (beta_axis * b_axis + 1 - beta_axis) * Eigen::Matrix3d::Identity();
  A_attract.bottomRightCorner(3, 3) = A_orientation;
  
  return  weight_attractor * A_attract;
}

//RMP for global damping
std::pair<Eigen::VectorXd, Eigen::MatrixXd> RiemannianMotionPolicy::calculate_global_damping(const Eigen::MatrixXd& Jp_obstacle) {
  Eigen::VectorXd f_damping = Eigen::VectorXd::Zero(6);
  Eigen::VectorXd velocity = Jp_obstacle * dq_;
  f_damping.topRows(3) = -k_damp * velocity * velocity.norm();
  Eigen::MatrixXd A_damping = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd identity_3 = Eigen::Matrix3d::Identity();
  A_damping.topLeftCorner(3, 3) = velocity.norm() * identity_3 * weight_damping;
  //return A_damping and f_damping
  return std::make_pair(f_damping, A_damping);
}

//RMP calculation of joint limit avoidance
void RiemannianMotionPolicy::rmp_joint_limit_avoidance(){
  //TODO: Implement the calculation of D_sigma fro joint limits
  //calculate sigma_u = 1/(1 + exp(-q))
  Eigen::VectorXd x_lower = Eigen::VectorXd::Zero(7);
  Eigen::VectorXd dx_lower = Eigen::VectorXd::Zero(7);
  Eigen::VectorXd x_upper = Eigen::VectorXd::Zero(7);
  Eigen::VectorXd dx_upper = Eigen::VectorXd::Zero(7);
  for (size_t i = 0; i < 7; ++i) {
    x_lower(i) = (q_(i) - q_lower_limit(i)) / (q_upper_limit(i) - q_lower_limit(i));
    dx_lower(i) = dq_(i) / (q_upper_limit(i) - q_lower_limit(i));

    x_upper(i) = (q_upper_limit(i) - q_(i)) / (q_upper_limit(i) - q_lower_limit(i));
    dx_upper(i) = -dq_(i) / (q_upper_limit(i) - q_lower_limit(i));

    f_joint_limits_lower(i) = kp_joint_limits/((std::pow(x_lower(i),2)/std::pow(l_p,2)) + accel_eps)  - kd_joint_limits * dx_lower(i);
    f_joint_limits_upper(i) = kp_joint_limits/((std::pow(x_upper(i),2)/std::pow(l_p,2)) + accel_eps)  - kd_joint_limits * dx_upper(i);
    A_joint_limits_lower(i,i) = weight_joint_limits * (1 - (1/(1 + exp(-dx_lower(i)/v_m)))) * (1/((x_lower(i)/l_m)+ epsilon_joint_limits));
    A_joint_limits_upper(i,i) = weight_joint_limits * (1 - (1/(1 + exp(-dx_upper(i)/v_m)))) * (1/((x_upper(i)/l_m)+ epsilon_joint_limits));
   
  }
}
 //RMP for joint limit velocity
 void RiemannianMotionPolicy::rmp_joint_velocity_limits(){
  Eigen::VectorXd q_dot_max = Eigen::VectorXd::Zero(7);
  q_dot_max << 2.175, 2.175, 2.175, 2.175, 2.61, 2.61, 2.61;
  
  for (size_t i = 0; i < 7; ++i) {
    double sign = std::copysign(1.0, dq_(i)); // Returns 1.0 for positive x, -1.0 for negative x.
    double dq_abs = std::abs(dq_(i));
    f_joint_velocity(i) = -k_joint_velocity * sign * (dq_abs - (q_dot_max(i) - 1.5));
    if(dq_abs < (q_dot_max(i) - 1.5)){
      A_joint_velocity(i,i) = 0.0;
    }
    else{
      A_joint_velocity(i,i) = weight_joint_velocity/(1 - std::pow((dq_abs - (q_dot_max(i) - 1.5)), 2)/(std::pow(1.5, 2)));
    }
  }
}

void RiemannianMotionPolicy::rmp_cspacetarget(){
  for (size_t i = 0; i < 7; ++i) {
    if(std::abs(q_(i)) < theta_cspace){
      f_c_space_target(i) = kp_c_space_target * (q_0(i) - q_(i)) - kd_c_space_target * dq_(i);
    }
    else{
      double diff = q_0(i) - q_(i);  // Scalar difference
      double abs_diff = std::abs(diff) + 1e-6;
      f_c_space_target(i) = kp_c_space_target * theta_cspace * diff/abs_diff - kd_c_space_target * dq_(i);
    }
  }
    A_c_space_target = weight_c_space_target * Eigen::MatrixXd::Identity(7,7);
}
//get global joint acceleration for torque calculation
void RiemannianMotionPolicy::get_ddq(){


  Eigen::MatrixXd I_77 = Eigen::MatrixXd::Identity(7, 7);
  Eigen::MatrixXd I_66 = Eigen::MatrixXd::Identity(6, 6);
  Eigen::MatrixXd Link2_a = jacobian2_obstacle.transpose() * A_obs_tilde2 * jacobian2_obstacle + jacobian2_obstacle.transpose() * A_damping2 * jacobian2_obstacle;
  Eigen::MatrixXd Link3_a = jacobian3_obstacle.transpose() * A_obs_tilde3 * jacobian3_obstacle + jacobian3_obstacle.transpose() * A_damping3 * jacobian3_obstacle;
  Eigen::MatrixXd Link4_a = jacobian4_obstacle.transpose() * A_obs_tilde4 * jacobian4_obstacle + jacobian4_obstacle.transpose() * A_damping4 * jacobian4_obstacle;
  Eigen::MatrixXd Link5_a = jacobian5_obstacle.transpose() * A_obs_tilde5 * jacobian5_obstacle + jacobian5_obstacle.transpose() * A_damping5 * jacobian5_obstacle;
  Eigen::MatrixXd Link6_a = jacobian6_obstacle.transpose() * A_obs_tilde6 * jacobian6_obstacle + jacobian6_obstacle.transpose() * A_damping6 * jacobian6_obstacle;
  Eigen::MatrixXd Link7_a = jacobian7_obstacle.transpose() * A_obs_tilde7 * jacobian7_obstacle + jacobian7_obstacle.transpose() * A_damping7 * jacobian7_obstacle;
  Eigen::MatrixXd Hand_a = jacobianhand_obstacle.transpose() * A_obs_tildehand * jacobianhand_obstacle + jacobianhand_obstacle.transpose() * A_dampinghand * jacobianhand_obstacle;
  Eigen::MatrixXd EE_a = jacobianEE_obstacle.transpose() * A_obs_tildeEE* jacobianEE_obstacle + jacobianEE_obstacle.transpose() * A_dampingEE * jacobianEE_obstacle;
  Eigen::MatrixXd Link2_b = jacobian2_obstacle.transpose() * A_obs_tilde2 * f_obs_tilde2 + jacobian2_obstacle.transpose() * A_damping2 * f_damping2;
  Eigen::MatrixXd Link3_b = jacobian3_obstacle.transpose() * A_obs_tilde3 * f_obs_tilde3 + jacobian3_obstacle.transpose() * A_damping3 * f_damping3; 
  Eigen::MatrixXd Link4_b = jacobian4_obstacle.transpose() * A_obs_tilde4 * f_obs_tilde4 + jacobian4_obstacle.transpose() * A_damping4 * f_damping4;
  Eigen::MatrixXd Link5_b = jacobian5_obstacle.transpose() * A_obs_tilde5 * f_obs_tilde5 + jacobian5_obstacle.transpose() * A_damping5 * f_damping5;
  Eigen::MatrixXd Link6_b = jacobian6_obstacle.transpose() * A_obs_tilde6 * f_obs_tilde6 + jacobian6_obstacle.transpose() * A_damping6 * f_damping6;
  Eigen::MatrixXd Link7_b = jacobian7_obstacle.transpose() * A_obs_tilde7 * f_obs_tilde7 + jacobian7_obstacle.transpose() * A_damping7 * f_damping7;
  Eigen::MatrixXd Hand_b = jacobianhand_obstacle.transpose() * A_obs_tildehand * f_obs_tildehand + jacobianhand_obstacle.transpose() * A_dampinghand * f_dampinghand;
  Eigen::MatrixXd EE_b = jacobianEE_obstacle.transpose() * A_obs_tildeEE * f_obs_tildeEE + jacobianEE_obstacle.transpose() * A_dampingEE * f_dampingEE;
  Eigen::MatrixXd A_total = jacobian.transpose()*A_attract *jacobian + Hand_a +EE_a +Link2_a + Link3_a + Link4_a + Link5_a + Link6_a + Link7_a + A_joint_limits_upper + A_joint_limits_lower + A_joint_velocity + A_c_space_target;
  Eigen::MatrixXd A_total_inv;
  pseudoInverse(A_total, A_total_inv); // get pseudoinverse for pullback 
  ddq_ =  A_total_inv* (jacobian.transpose() * A_attract * x_dd_des + Hand_b +EE_b +Link2_b + Link3_b + Link4_b + Link5_b + Link6_b + Link7_b 
                          + A_joint_limits_upper * f_joint_limits_upper + A_joint_limits_lower * f_joint_limits_lower
                          + A_joint_velocity * f_joint_velocity + A_c_space_target * f_c_space_target);
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
  //Load parameters from yaml file
  YAML::Node config = YAML::LoadFile("/home/andri/franka_ros2_ws/src/riemannian_motion_policy/src/config.yaml");

  // Load obstacle avoidance parameters
  eta_rep = config["obstacle_avoidance"]["eta_rep"].as<double>();
  mu_rep = config["obstacle_avoidance"]["mu_rep"].as<double>();
  eta_damp = config["obstacle_avoidance"]["eta_damp"].as<double>();
  mu_damp = config["obstacle_avoidance"]["mu_damp"].as<double>();
  epsilon = config["obstacle_avoidance"]["epsilon"].as<double>();
  weight_obstacle = config["obstacle_avoidance"]["weight_obstacle"].as<double>();

  // Load attractor parameters
  alpha_min = config["attractor"]["alpha_min"].as<double>();
  sigma_a = config["attractor"]["sigma_a"].as<double>();
  sigma_b = config["attractor"]["sigma_b"].as<double>();
  b = config["attractor"]["b"].as<double>();
  sigma_o = config["attractor"]["sigma_o"].as<double>();
  b_axis = config["attractor"]["b_axis"].as<double>();
  weight_attractor = config["attractor"]["weight_attractor"].as<double>();

  // Load global damping parameters
  k_damp = config["global_damping"]["k_damp"].as<double>();
  weight_damping = config["global_damping"]["weight_damping"].as<double>();

  // Load velocity limits parameters
  k_joint_velocity = config["velocity_limits"]["k_joint_velocity"].as<double>();
  weight_joint_velocity = config["velocity_limits"]["weight_joint_velocity"].as<double>();

  // Load joint limits parameters
  kp_joint_limits = config["joint_limits"]["kp_joint_limits"].as<double>();
  kd_joint_limits = config["joint_limits"]["kd_joint_limits"].as<double>();
  l_m = config["joint_limits"]["metric_length_scale"].as<double>();
  epsilon_joint_limits = config["joint_limits"]["epsilon_joint_limits"].as<double>();
  v_m = config["joint_limits"]["metric_velocity_length_scale"].as<double>();
  l_p = config["joint_limits"]["accel_exploder_length_scale"].as<double>();
  accel_eps = config["joint_limits"]["accel_eps"].as<double>();
  weight_joint_limits = config["joint_limits"]["weight_joint_limits"].as<double>();

  // Load C-space Target parameters
  kp_c_space_target = config["c_space_target"]["kp_c_space_target"].as<double>();
  kd_c_space_target = config["c_space_target"]["kd_c_space_target"].as<double>();
  theta_cspace = config["c_space_target"]["theta"].as<double>();
  weight_c_space_target = config["c_space_target"]["weight_c_space_target"].as<double>();


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
    d_obshand << msg->framehandx, msg->framehandy, msg->framehandz;
    d_obsEE << msg->frameeex, msg->frameeey, msg->frameeez;
    // Handle the Jacobian of the closest point
    jacobian_array2 = msg->jacobian2;
    jacobian_array3 = msg->jacobian3;
    jacobian_array4 = msg->jacobian4;
    jacobian_array5 = msg->jacobian5;
    jacobian_array6 = msg->jacobian6;
    jacobian_array7 = msg->jacobian7;
    jacobian_arrayhand = msg->jacobianhand;
    jacobian_arrayEE = msg->jacobianee;
    //reshape to matrix (6x7) column major
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::ColMajor>> jacobian2obstacle(jacobian_array2.data());
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::ColMajor>> jacobian3obstacle(jacobian_array3.data());
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::ColMajor>> jacobian4obstacle(jacobian_array4.data());
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::ColMajor>> jacobian5obstacle(jacobian_array5.data());
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::ColMajor>> jacobian6obstacle(jacobian_array6.data());
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::ColMajor>> jacobian7obstacle(jacobian_array7.data());
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::ColMajor>> jacobianhandobstacle(jacobian_arrayhand.data());
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::ColMajor>> jacobianEEobstacle(jacobian_arrayEE.data());
    //assign to class variables
    jacobian2_obstacle = jacobian2obstacle;
    jacobian3_obstacle = jacobian3obstacle;
    jacobian4_obstacle = jacobian4obstacle;
    jacobian5_obstacle = jacobian5obstacle;
    jacobian6_obstacle = jacobian6obstacle;
    jacobian7_obstacle = jacobian7obstacle;
    jacobianhand_obstacle = jacobianhandobstacle;
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
  
  Jp_obstacle2 = jacobian2_obstacle.topRows(3);
  Jp_obstacle3 = jacobian3_obstacle.topRows(3);
  Jp_obstacle4 = jacobian4_obstacle.topRows(3);
  Jp_obstacle5 = jacobian5_obstacle.topRows(3);
  Jp_obstacle6 = jacobian6_obstacle.topRows(3);
  Jp_obstacle7 = jacobian7_obstacle.topRows(3);
  Jp_obstaclehand = jacobianhand_obstacle.topRows(3);
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


  //d_obs1 = calculateNearestPointOnSphere(position, sphere_center, sphere_radius);
  //d_obs1 = d_obs_prev1 * 0.99 + d_obs1 * 0.01;
  Lambda = (jacobian * M.inverse() * jacobian.transpose()).inverse();
  x_dd_des = (-K_RMP * (error) - D_RMP* jacobian * dq_);
  A_attract = calculate_target_attraction(error, jacobian);
  //A_attract = Eigen::MatrixXd::Zero(6, 6);
  f_obs_tildeEE = calculate_f_obstacle(d_obsEE, Jp_obstacleEE);
  A_obs_tildeEE = calculate_A_obstacle(d_obsEE, f_obs_tildeEE, 0.5, Jp_obstacleEE);
  f_obs_tilde2 = calculate_f_obstacle(d_obs2, Jp_obstacle2);
  A_obs_tilde2 = calculate_A_obstacle(d_obs2, f_obs_tilde2, 0.5, Jp_obstacle2);
  f_obs_tilde3 = calculate_f_obstacle(d_obs3, Jp_obstacle3);
  A_obs_tilde3 = calculate_A_obstacle(d_obs3, f_obs_tilde3, 0.5, Jp_obstacle3);
  f_obs_tilde4 = calculate_f_obstacle(d_obs4, Jp_obstacle4);
  A_obs_tilde4 = calculate_A_obstacle(d_obs4, f_obs_tilde4, 0.5,  Jp_obstacle4);
  f_obs_tilde5 = calculate_f_obstacle(d_obs5, Jp_obstacle5);
  A_obs_tilde5 = calculate_A_obstacle(d_obs5, f_obs_tilde5, 0.5, Jp_obstacle5);
  f_obs_tilde6 = calculate_f_obstacle(d_obs6, Jp_obstacle6);
  A_obs_tilde6 = calculate_A_obstacle(d_obs6, f_obs_tilde6, 0.5,   Jp_obstacle6);
  f_obs_tilde7 = calculate_f_obstacle(d_obs7, Jp_obstacle7);
  A_obs_tilde7 = calculate_A_obstacle(d_obs7, f_obs_tilde7, 0.5, Jp_obstacle7);
  f_obs_tildehand = calculate_f_obstacle(d_obshand, Jp_obstaclehand);
  A_obs_tildehand = calculate_A_obstacle(d_obshand, f_obs_tildehand, 0.5, Jp_obstaclehand);
  auto [f_damping2, A_damping2] = calculate_global_damping(Jp_obstacle2);
  auto [f_damping3, A_damping3] = calculate_global_damping(Jp_obstacle3);
  auto [f_damping4, A_damping4] = calculate_global_damping(Jp_obstacle4);
  auto [f_damping5, A_damping5] = calculate_global_damping(Jp_obstacle5);
  auto [f_damping6, A_damping6] = calculate_global_damping(Jp_obstacle6);
  auto [f_damping7, A_damping7] = calculate_global_damping(Jp_obstacle7);
  auto [f_dampinghand, A_dampinghand] = calculate_global_damping(Jp_obstaclehand);
  auto [f_dampingEE, A_dampingEE] = calculate_global_damping(Jp_obstacleEE);
  rmp_joint_limit_avoidance();
  rmp_joint_velocity_limits();
  rmp_cspacetarget();
  get_ddq();
  
  // Calculate the desired torque
  tau_RMP = M * ddq_;
  // Calculate friction torques
  //calculate_tau_friction();
  calculate_tau_gravity(coriolis, gravity_force_vector, jacobian);
  //tau_gravity_error = tau_gravity - gravity_force_vector;

  auto tau_d_placeholder = tau_RMP + coriolis; //add nullspace, friction, gravity and coriolis components to desired torque
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
    std::cout<<"ddq" << std::endl;
    std::cout << ddq_<< std::endl;
    std::cout<<"f_obs_tildeEE" << std::endl;
    std::cout << f_obs_tildeEE << std::endl;
    std::cout<<"f_joint_velocity" << std::endl;
    std::cout << f_joint_velocity << std::endl;
    std::cout<<"f_dampingEE" << std::endl;
    std::cout << f_dampingEE << std::endl;
    //std::cout << "error_pose" << std::endl;
    //std::cout << error << std::endl;
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