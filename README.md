# Riemannian Motion Policy
### ROS2 RMP_controller

Prerequisites:
* ROS2 humble <br />
* Libfranka 0.13.0 or newer <br />
* franka_ros2 v0.13.1 <br />

For further information, please refer to the [Franka ROS2 FCI documentation](https://support.franka.de/docs/franka_ros2.html)

Once you have everything set up, follow the steps below to get the controller running.

### RMP installation
Clone this repository in the src directory of your franka_ros2_ws: <br />
```bash
cd franka_ros2_ws/src 
git clone https://github.com/acaviezel/Riemannian-Motion-Policies-Franka-Emika-Robot.git
```
For the moment, you need to add the following lines of code, to your controllers.yaml file inside franka_ros2/franka_bringup/config/:
```bash
riemannian_motion_policy:
      type: riemannian_motion_policy/RiemannianMotionPolicy
```

Clone the messages package in the src directory: <br />
```bash
cd franka_ros2_ws/src
git clone https://github.com/acaviezel/messages_fr3.git
```

Build the package or whole workspace: <br />
```bash
colcon build --packages-select cartesian_impedance_control --cmake-args -DCMAKE_BUILD_TYPE=Release
colcon build --cmake-args -DCMAKE_BUILD_TYPE=Release #Builds all the packages in your src folder
```
If not yet done, ensure your setup is always source by adding the following line to the end of your .bashrc file (to get access to it, you need to execute `nano .bashrc` in your home directory). : <br />
```bash
source /home/<user>/franka_ros2_ws/install/setup.sh 
```

### Distance Calculator and MoveIt Scene Loader
Follow the installation guidelines from [motion_planning_mt](https://github.com/acaviezel/motion_planning_mt)

### Launch controller with MoveIt
After building the whole workspace, the controller can now be launched.
1. Launch the moveit environment and the RMP controller.
```
ros2 launch franka_moveit_config moveit.launch.py robot_ip:=<fci-ip>
```
2. Run the scene node, which populates three cylinders (adjustable).
```
ros2 run motion_planning_mt cylinder_scene
```
3. Run the distance calculator node, which calculates the distance between the closest obstacle and each robot link.
```
ros2 run motion_planning_mt distance_calculator
```
To visualize the obstacles and the minimum distance in RViz, add two ROS2 topics. In displays panel, click "Add" button on the lower left corner.
   1. Add one MarkerArray, set the topic to which it subscribes to `/rviz_visual_tools` (If not already subscribed).
   2. Add one MarkerArray, set the topic to which it subscribes to `/minimum_distance_visualization`.

