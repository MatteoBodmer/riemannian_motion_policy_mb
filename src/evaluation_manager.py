
import rclpy
from rclpy.node import Node
from geometry_msgs.msg import Pose
from moveit_msgs.msg import CollisionObject, PlanningScene
from shape_msgs.msg import SolidPrimitive
from std_msgs.msg import Header
import random
import time
import json
import numpy as np
from datetime import datetime
import threading

class EvaluationManager(Node):
    def __init__(self):
        super().__init__('evaluation_manager')
        
        # Publishers
        self.pose_pub = self.create_publisher(Pose, '/riemannian_motion_policy/reference_pose', 10)
        self.planning_scene_pub = self.create_publisher(PlanningScene, '/planning_scene', 10)
        
        # Evaluation parameters
        self.num_simulations = 5  # Start with smaller number for testing
        self.current_simulation = 0
        self.results = []
        
        # Workspace bounds (adjust these to your robot's workspace)
        self.workspace_bounds = {
            'x': (0.2, 0.8),
            'y': (-0.5, 0.5), 
            'z': (0.2, 0.8)
        }
        
        # Wait for other nodes to be ready
        time.sleep(5.0)
        
        self.get_logger().info("Evaluation Manager initialized. Starting automated evaluation...")
        
        # Start evaluation in a separate thread
        self.eval_thread = threading.Thread(target=self.run_evaluation_suite)
        self.eval_thread.start()
    
    def generate_random_obstacle_pose(self):
        """Generate random obstacle position within workspace bounds"""
        pose = Pose()
        pose.position.x = random.uniform(*self.workspace_bounds['x'])
        pose.position.y = random.uniform(*self.workspace_bounds['y'])
        pose.position.z = random.uniform(*self.workspace_bounds['z'])
        
        # Keep orientation upright for simplicity
        pose.orientation.w = 1.0
        pose.orientation.x = 0.0
        pose.orientation.y = 0.0
        pose.orientation.z = 0.0
        
        return pose
    
    def generate_random_target_pose(self, obstacle_pose):
        """Generate random target pose avoiding obstacle vicinity"""
        min_distance = 0.15  # 15cm clearance from obstacle
        max_attempts = 50
        
        for _ in range(max_attempts):
            target = Pose()
            target.position.x = random.uniform(*self.workspace_bounds['x'])
            target.position.y = random.uniform(*self.workspace_bounds['y'])
            target.position.z = random.uniform(*self.workspace_bounds['z'])
            
            # Keep orientation upright
            target.orientation.w = 1.0
            target.orientation.x = 0.0
            target.orientation.y = 0.0
            target.orientation.z = 0.0
            
            # Check distance from obstacle
            distance = np.sqrt(
                (target.position.x - obstacle_pose.position.x)**2 +
                (target.position.y - obstacle_pose.position.y)**2 +
                (target.position.z - obstacle_pose.position.z)**2
            )
            
            if distance > min_distance:
                return target
        
        # If no valid target found, use a safe default
        target = Pose()
        target.position.x = 0.3
        target.position.y = 0.0
        target.position.z = 0.6
        target.orientation.w = 1.0
        return target
    
    def spawn_obstacle(self, obstacle_pose):
        """Spawn a cylinder obstacle at given pose"""
        # Clear existing obstacles first
        self.clear_obstacles()
        
        # Create collision object
        collision_obj = CollisionObject()
        collision_obj.header.frame_id = "fr3_link0"
        collision_obj.header.stamp = self.get_clock().now().to_msg()
        collision_obj.id = "evaluation_obstacle"
        
        # Create cylinder primitive
        cylinder = SolidPrimitive()
        cylinder.type = SolidPrimitive.CYLINDER
        cylinder.dimensions = [0.4, 0.05]  # height=0.4m, radius=0.05m
        
        collision_obj.primitives.append(cylinder)
        collision_obj.primitive_poses.append(obstacle_pose)
        collision_obj.operation = CollisionObject.ADD
        
        # Create planning scene message
        planning_scene = PlanningScene()
        planning_scene.world.collision_objects.append(collision_obj)
        planning_scene.is_diff = True
        
        # Publish multiple times to ensure reception
        for _ in range(5):
            self.planning_scene_pub.publish(planning_scene)
            time.sleep(0.1)
        
        self.get_logger().info(f"Spawned obstacle at ({obstacle_pose.position.x:.2f}, {obstacle_pose.position.y:.2f}, {obstacle_pose.position.z:.2f})")
    
    def clear_obstacles(self):
        """Clear all existing obstacles"""
        collision_obj = CollisionObject()
        collision_obj.header.frame_id = "fr3_link0"
        collision_obj.header.stamp = self.get_clock().now().to_msg()
        collision_obj.id = "evaluation_obstacle"
        collision_obj.operation = CollisionObject.REMOVE
        
        planning_scene = PlanningScene()
        planning_scene.world.collision_objects.append(collision_obj)
        planning_scene.is_diff = True
        
        for _ in range(3):
            self.planning_scene_pub.publish(planning_scene)
            time.sleep(0.1)
    
    def set_target_pose(self, target_pose):
        """Send target pose to RMP controller"""
        for _ in range(5):  # Publish multiple times to ensure reception
            self.pose_pub.publish(target_pose)
            time.sleep(0.1)
        
        self.get_logger().info(f"Set target pose: ({target_pose.position.x:.2f}, {target_pose.position.y:.2f}, {target_pose.position.z:.2f})")
    
    def run_single_simulation(self):
        """Run a single simulation trial"""
        self.get_logger().info(f"Starting simulation {self.current_simulation + 1}/{self.num_simulations}")
        
        # 1. Generate random obstacle
        obstacle_pose = self.generate_random_obstacle_pose()
        self.spawn_obstacle(obstacle_pose)
        
        # 2. Generate random target pose
        target_pose = self.generate_random_target_pose(obstacle_pose)
        
        # 3. Wait for obstacle to be registered
        time.sleep(2.0)
        
        # 4. Set target pose
        self.set_target_pose(target_pose)
        
        # 5. Monitor execution for fixed time
        execution_time = 10.0  # seconds
        self.get_logger().info(f"Monitoring execution for {execution_time} seconds...")
        time.sleep(execution_time)
        
        # 6. Record results (basic for now)
        result = {
            'simulation_id': self.current_simulation,
            'obstacle_position': [obstacle_pose.position.x, obstacle_pose.position.y, obstacle_pose.position.z],
            'target_position': [target_pose.position.x, target_pose.position.y, target_pose.position.z],
            'execution_time': execution_time,
            'timestamp': datetime.now().isoformat()
        }
        
        return result
    
    def reset_simulation(self):
        """Reset simulation for next trial"""
        # Clear obstacles
        self.clear_obstacles()
        
        # Reset to home position (you might want to add a service for this)
        home_pose = Pose()
        home_pose.position.x = 0.3
        home_pose.position.y = 0.0
        home_pose.position.z = 0.6
        home_pose.orientation.w = 1.0
        self.set_target_pose(home_pose)
        
        # Wait for reset
        time.sleep(5.0)
    
    def run_evaluation_suite(self):
        """Run the complete evaluation suite"""
        self.get_logger().info(f"Starting evaluation suite with {self.num_simulations} simulations")
        
        for i in range(self.num_simulations):
            self.current_simulation = i
            
            try:
                result = self.run_single_simulation()
                self.results.append(result)
                self.get_logger().info(f"Completed simulation {i+1}/{self.num_simulations}")
                
                # Reset for next simulation (except for the last one)
                if i < self.num_simulations - 1:
                    self.reset_simulation()
                    
            except Exception as e:
                self.get_logger().error(f"Error in simulation {i+1}: {str(e)}")
                self.reset_simulation()
        
        # Save results
        self.save_results()
        self.get_logger().info("Evaluation suite completed!")
    
    def save_results(self):
        """Save evaluation results to file"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"/home/matteo/franka_ros2_ws/evaluation_results_{timestamp}.json"
        
        summary = {
            'total_simulations': len(self.results),
            'timestamp': datetime.now().isoformat(),
            'results': self.results
        }
        
        with open(filename, 'w') as f:
            json.dump(summary, f, indent=2)
        
        self.get_logger().info(f"Results saved to {filename}")

def main(args=None):
    rclpy.init(args=args)
    
    evaluation_manager = EvaluationManager()
    
    try:
        rclpy.spin(evaluation_manager)
    except KeyboardInterrupt:
        pass
    finally:
        evaluation_manager.destroy_node()
        rclpy.shutdown()

if __name__ == '__main__':
    main()