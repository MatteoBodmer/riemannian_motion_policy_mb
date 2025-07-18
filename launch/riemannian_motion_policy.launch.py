import os
import xacro

from ament_index_python.packages import get_package_share_directory

from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument, OpaqueFunction, ExecuteProcess, RegisterEventHandler
from launch.event_handlers import OnProcessExit

from launch.actions import IncludeLaunchDescription
from launch.launch_description_sources import PythonLaunchDescriptionSource
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node


def get_robot_description(context, arm_id, load_gripper, franka_hand):
    arm_id_str = context.perform_substitution(arm_id)
    load_gripper_str = context.perform_substitution(load_gripper)
    franka_hand_str = context.perform_substitution(franka_hand)

    franka_xacro_file = os.path.join(
        get_package_share_directory('franka_description'),
        'robots',
        arm_id_str,
        arm_id_str + '.urdf.xacro'
    )

    robot_description_config = xacro.process_file(
        franka_xacro_file,
        mappings={
            'arm_id': arm_id_str,
            'hand': load_gripper_str,
            'ros2_control': 'true',
            'gazebo': 'true',
            'ee_id': franka_hand_str,
            'gazebo_effort': 'true'
        }
    )
    robot_description = {'robot_description': robot_description_config.toxml()}

    robot_state_publisher = Node(
        package='robot_state_publisher',
        executable='robot_state_publisher',
        name='robot_state_publisher',
        output='screen',
        parameters=[
            robot_description,
        ]
    )

    return [robot_state_publisher]


def generate_launch_description():
    # Configure ROS nodes for launch
    load_gripper_name = 'load_gripper'
    franka_hand_name = 'franka_hand'
    arm_id_name = 'arm_id'

    load_gripper = LaunchConfiguration(load_gripper_name)
    franka_hand = LaunchConfiguration(franka_hand_name)
    arm_id = LaunchConfiguration(arm_id_name)

    load_gripper_launch_argument = DeclareLaunchArgument(
        load_gripper_name,
        default_value='true',
        description='true/false for activating the gripper')
    franka_hand_launch_argument = DeclareLaunchArgument(
        franka_hand_name,
        default_value='franka_hand',
        description='Default value: franka_hand')
    arm_id_launch_argument = DeclareLaunchArgument(
        arm_id_name,
        default_value='fr3',
        description='Available values: fr3, fp3 and fer')

    # Get robot description
    robot_state_publisher = OpaqueFunction(
        function=get_robot_description,
        args=[arm_id, load_gripper, franka_hand])

    # Gazebo Sim
    os.environ['GZ_SIM_RESOURCE_PATH'] = os.path.dirname(get_package_share_directory('franka_description'))
    pkg_ros_gz_sim = get_package_share_directory('ros_gz_sim')
    gazebo_empty_world = IncludeLaunchDescription(
        PythonLaunchDescriptionSource(
            os.path.join(pkg_ros_gz_sim, 'launch', 'gz_sim.launch.py')),
        launch_arguments={'gz_args': 'empty.sdf -r', }.items(),
    )

    # Spawn the Franka robot in Gazebo
    spawn_robot = Node(
        package='ros_gz_sim',
        executable='create',
        arguments=['-topic', '/robot_description'],
        output='screen',
    )

    # Define the absolute path to the obstacle URDF
    obstacle_urdf_path = '/home/matteo/franka_ros2_ws/src/Riemannian-Motion-Policies-Franka-Emika-Robot/src/Objects.urdf'

    spawn_obstacle = Node(
        package='ros_gz_sim',
        executable='create',
        arguments=[
            '-file', obstacle_urdf_path,
            '-name', 'obstacle',
            '-x', '0.5', '-y', '0.0', '-z', '0.5'
        ],
        output='screen',
    )


    # Visualize in RViz
    # rviz_file = os.path.join(get_package_share_directory('franka_description'), 'rviz', 'visualize_franka.rviz')
    # rviz = Node(
    #     package='rviz2',
    #     executable='rviz2',
    #     name='rviz2',
    #     arguments=['--display-config', rviz_file, '-f', 'world'],
    # )

    load_joint_state_broadcaster = ExecuteProcess(
        cmd=['ros2', 'control', 'load_controller', '--set-state', 'active', 'joint_state_broadcaster'],
        output='screen'
    )

    riemannian_motion_policy = ExecuteProcess(
        cmd=['ros2', 'control', 'load_controller', '--set-state', 'active', 'riemannian_motion_policy'],
        output='screen'
    )

    set_load = ExecuteProcess(
        cmd=['/home/matteo/franka_ros2_ws/src/Riemannian-Motion-Policies-Franka-Emika-Robot/launch/set_load.sh'],
        output='screen',
    )

    start_controller = Node(
        package='controller_manager',
        executable='spawner',
        arguments=['riemannian_motion_policy'],
        output='screen',
    )


    return LaunchDescription([
        load_gripper_launch_argument,
        franka_hand_launch_argument,
        arm_id_launch_argument,
        gazebo_empty_world,
        robot_state_publisher,
        #rviz,
        spawn_robot,
        

        # spawn the obstacle after spawning the robot
        RegisterEventHandler(
            event_handler=OnProcessExit(
                target_action=spawn_robot,
                on_exit=[spawn_obstacle],
            )
        ),

        RegisterEventHandler(
            event_handler=OnProcessExit(
                target_action=spawn_obstacle,
                on_exit=[load_joint_state_broadcaster],
            )
        ),

        RegisterEventHandler(
            event_handler=OnProcessExit(
                target_action=load_joint_state_broadcaster,
                on_exit=[riemannian_motion_policy],
            )
        ),

        Node(
            package='joint_state_publisher',
            executable='joint_state_publisher',
            name='joint_state_publisher',
            parameters=[
                {'source_list': ['joint_states'],
                 'rate': 30}],
        ),

        set_load,

        RegisterEventHandler(
            event_handler=OnProcessExit(
                target_action=set_load,
                on_exit=[start_controller],
            )
        ),
    ])
