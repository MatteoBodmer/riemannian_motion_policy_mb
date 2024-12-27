# Riemannian Motion Policy (RMP) Controller

### Prerequisites
Ensure you have the following installed and configured:

- **ROS 2 Humble**
- **libfranka v0.13.0 or newer**
- **franka_ros2 v0.13.1 or newer**

For detailed setup instructions, please refer to the [Franka ROS2 FCI documentation](https://support.franka.de/docs/franka_ros2.html).

---

## Installation

### Step 1: Clone the RMP Repository
Clone this repository into the `src` directory of your `franka_ros2_ws` workspace:

```bash
cd ~/franka_ros2_ws/src
git clone https://github.com/acaviezel/Riemannian-Motion-Policies-Franka-Emika-Robot.git

Step 2: Modify Configuration for RMP Controller

Add the following lines to your controllers.yaml file, located at franka_ros2/franka_bringup/config/:

riemannian_motion_policy:
  type: riemannian_motion_policy/RiemannianMotionPolicy

Step 3: Clone the Messages Package

Clone the messages_fr3 package into the src directory of your workspace:

cd ~/franka_ros2_ws/src
git clone https://github.com/acaviezel/messages_fr3.git

Step 4: Build the Workspace

Build either the specific package or the entire workspace:

    To build the RMP controller package:

colcon build --packages-select cartesian_impedance_control --cmake-args -DCMAKE_BUILD_TYPE=Release

To build all packages in the src folder:

    colcon build --cmake-args -DCMAKE_BUILD_TYPE=Release

Step 5: Update .bashrc File

To ensure your setup is sourced each time you open a terminal, add the following line to the end of your .bashrc file:

source ~/franka_ros2_ws/install/setup.sh

To edit the .bashrc file:

nano ~/.bashrc

After editing, reload the file:

source ~/.bashrc

Additional Setup: Distance Calculator and MoveIt Scene Loader

Follow the installation instructions provided in the motion_planning_mt repository.
Launching the RMP Controller with MoveIt
Step 1: Start the MoveIt Environment

Launch the MoveIt environment and the RMP controller:

ros2 launch franka_moveit_config moveit.launch.py robot_ip:=<fci-ip>

Replace <fci-ip> with the IP address of your Franka Control Interface (FCI).
Step 2: Launch the Scene Node

Run the scene node to populate the environment with three cylinders (adjustable as needed):

ros2 run motion_planning_mt cylinder_scene

Step 3: Launch the Distance Calculator

Start the node to calculate the distance between each robot link and the nearest obstacle:

ros2 run motion_planning_mt distance_calculator

Step 4: Visualize Obstacles and Minimum Distances in RViz

To visualize the obstacles and the minimum distance in RViz:

    Open the Displays panel in RViz.
    Click the "Add" button (bottom left corner).
    Add the following:
        MarkerArray: Set the topic to /rviz_visual_tools.
        MarkerArray: Set the topic to /minimum_distance_visualization.

