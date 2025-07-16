import os
import xacro
import yaml

from ament_index_python.packages import get_package_share_directory

from launch import LaunchDescription
from launch.actions import (
    DeclareLaunchArgument, ExecuteProcess, OpaqueFunction,
    IncludeLaunchDescription, RegisterEventHandler
)
from launch.event_handlers import OnProcessExit
from launch.launch_description_sources import PythonLaunchDescriptionSource
from launch.substitutions import LaunchConfiguration, Command, FindExecutable
from launch_ros.actions import Node
from launch_ros.parameter_descriptions import ParameterValue


def load_yaml(package_name, file_path):
    package_path = get_package_share_directory(package_name)
    absolute_file_path = os.path.join(package_path, file_path)

    try:
        with open(absolute_file_path, 'r') as file:
            return yaml.safe_load(file)
    except EnvironmentError:
        return None


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
        parameters=[robot_description],
        remappings=[('joint_states', '/gazebo/joint_states')],
    )

    return [robot_state_publisher]


def generate_launch_description():
    # Common arguments
    load_gripper_arg = DeclareLaunchArgument('load_gripper', default_value='true')
    franka_hand_arg = DeclareLaunchArgument('franka_hand', default_value='franka_hand')
    arm_id_arg = DeclareLaunchArgument('arm_id', default_value='fr3')
    robot_ip_arg = DeclareLaunchArgument('robot_ip', default_value='sim', description='Hostname or IP of robot.')
    use_fake_hw_arg = DeclareLaunchArgument('use_fake_hardware', default_value='true')
    fake_sensor_cmd_arg = DeclareLaunchArgument('fake_sensor_commands', default_value='true')
    db_arg = DeclareLaunchArgument('db', default_value='False')

    # Launch configurations
    load_gripper = LaunchConfiguration('load_gripper')
    franka_hand = LaunchConfiguration('franka_hand')
    arm_id = LaunchConfiguration('arm_id')
    robot_ip = LaunchConfiguration('robot_ip')
    use_fake_hw = LaunchConfiguration('use_fake_hardware')
    fake_sensor_cmd = LaunchConfiguration('fake_sensor_commands')

    # === RMP & Gazebo ===
    os.environ['GZ_SIM_RESOURCE_PATH'] = os.path.dirname(get_package_share_directory('franka_description'))
    pkg_ros_gz_sim = get_package_share_directory('ros_gz_sim')

    gazebo_world = IncludeLaunchDescription(
        PythonLaunchDescriptionSource(
            os.path.join(pkg_ros_gz_sim, 'launch', 'gz_sim.launch.py')),
        launch_arguments={'gz_args': 'empty.sdf -r'}.items(),
    )

    robot_state_publisher = OpaqueFunction(
        function=get_robot_description,
        args=[arm_id, load_gripper, franka_hand]
    )

    spawn_robot = Node(
        package='ros_gz_sim',
        executable='create',
        arguments=['-topic', '/robot_description'],
        output='screen',
    )

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

    # === MoveIt & RViz ===

    franka_semantic_xacro_file = os.path.join(
        get_package_share_directory('franka_fr3_moveit_config'),
        'srdf', 'fr3_arm.srdf.xacro')

    robot_description_semantic_str = Command([
        FindExecutable(name='xacro'), ' ', franka_semantic_xacro_file, ' hand:=true'
    ])

    robot_description_semantic = {
        'robot_description_semantic': ParameterValue(
            robot_description_semantic_str,
            value_type=str)
    }

    kinematics_yaml = load_yaml('franka_fr3_moveit_config', 'config/kinematics.yaml')
    ompl_planning_yaml = load_yaml('franka_fr3_moveit_config', 'config/ompl_planning.yaml')
    ompl_planning_pipeline_config = {
        'move_group': {
            'planning_plugin': 'ompl_interface/OMPLPlanner',
            'request_adapters': 'default_planner_request_adapters/AddTimeOptimalParameterization '
                                'default_planner_request_adapters/ResolveConstraintFrames '
                                'default_planner_request_adapters/FixWorkspaceBounds '
                                'default_planner_request_adapters/FixStartStateBounds '
                                'default_planner_request_adapters/FixStartStateCollision '
                                'default_planner_request_adapters/FixStartStatePathConstraints',
            'start_state_max_bounds_error': 0.1,
        }
    }
    ompl_planning_pipeline_config['move_group'].update(ompl_planning_yaml)

    trajectory_execution = {
        'moveit_manage_controllers': True,
        'trajectory_execution.allowed_execution_duration_scaling': 1.2,
        'trajectory_execution.allowed_goal_duration_margin': 0.5,
        'trajectory_execution.allowed_start_tolerance': 0.01,
    }

    moveit_controllers_yaml = load_yaml('franka_fr3_moveit_config', 'config/fr3_controllers.yaml')
    moveit_controllers = {
        'moveit_simple_controller_manager': moveit_controllers_yaml,
        'moveit_controller_manager': 'moveit_simple_controller_manager/MoveItSimpleControllerManager',
    }

    planning_scene_monitor_parameters = {
        'publish_planning_scene': True,
        'publish_geometry_updates': True,
        'publish_state_updates': True,
        'publish_transforms_updates': True,
        'publish_robot_description_semantic' : True,
        'monitor_octomap': False,
    }

    # Extra node to publish robot_description_semantic globally
    semantic_param_publisher = Node(
        package='robot_state_publisher',
        executable='robot_state_publisher',
        name='semantic_param_publisher',
        output='screen',
        parameters=[robot_description_semantic],
    )

    run_move_group_node = Node(
        package='moveit_ros_move_group',
        executable='move_group',
        output='screen',
        parameters=[
            robot_description_semantic,
            kinematics_yaml,
            ompl_planning_pipeline_config,
            trajectory_execution,
            moveit_controllers,
            planning_scene_monitor_parameters,
        ],
    )

    rviz_config = os.path.join(
        get_package_share_directory('franka_fr3_moveit_config'),
        'rviz', 'moveit.rviz')

    rviz_node = Node(
        package='rviz2',
        executable='rviz2',
        name='rviz2',
        output='log',
        arguments=['-d', rviz_config],
        parameters=[
            robot_description_semantic,
            kinematics_yaml,
            ompl_planning_pipeline_config,
        ],
    )

    # === Final LaunchDescription ===

    return LaunchDescription([
        load_gripper_arg, franka_hand_arg, arm_id_arg,
        robot_ip_arg, use_fake_hw_arg, fake_sensor_cmd_arg, db_arg,

        gazebo_world,
        robot_state_publisher,
        spawn_robot,

        RegisterEventHandler(
            OnProcessExit(
                target_action=spawn_robot,
                on_exit=[spawn_obstacle],
            )
        ),
        RegisterEventHandler(
            OnProcessExit(
                target_action=spawn_obstacle,
                on_exit=[load_joint_state_broadcaster],
            )
        ),
        RegisterEventHandler(
            OnProcessExit(
                target_action=load_joint_state_broadcaster,
                on_exit=[riemannian_motion_policy],
            )
        ),
        set_load,
        RegisterEventHandler(
            OnProcessExit(
                target_action=set_load,
                on_exit=[start_controller],
            )
        ),

        # Publish robot_description_semantic globally
        semantic_param_publisher,

        # MoveIt & RViz
        run_move_group_node,
        rviz_node
    ])
