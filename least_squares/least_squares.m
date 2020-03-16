close all
#clear
#clc

#load dependencies
addpath "../"
addpath "../tools/g2o_wrapper"
addpath "../tools/visualization"
addpath "../ekf_slam"
source "../tools/utilities/geometry_helpers_2d.m"
source "./do_ls.m"

[_, poses, transitions, observations] = loadG2o("../datasets/BearingOnlySLAM/slam2D_bearing_only_initial_guess.g2o");
[lm_real, poses_real, transitions, observations] = loadG2o("../datasets/BearingOnlySLAM/slam2D_bearing_only_ground_truth.g2o");
#transitions and observations seem to be equal, poses differ

#initalization

[x_lm, id_to_guess, guess_to_id] = init(poses, observations);
#[x_lm, id_to_guess, guess_to_id] = init(poses_real, observations);
number_of_lm = rows(x_lm);

#draw true landmarks
if(false)
  figure(1)
  drawLandmarks(lm_real,'g', 'fill')
  lm_guess = landmark(guess_to_id(1), x_lm(1,:));
  for u = 2:rows(x_lm);
    lm_guess(end+1) = landmark(guess_to_id(u), x_lm(u,:));
  endfor
  drawLandmarks(lm_guess,'r', 'fill')

  #draw true trajectory
  figure(2)
  for u = 1:length(poses_real)-1
    #draw line from begin to end
    x_begin = poses_real(u).x;
    y_begin = poses_real(u).y;
    x_end   = poses_real(u+1).x;
    y_end   = poses_real(u+1).y;
    plot([x_begin x_end], [y_begin y_end], "g", "linewidth", 2);
    hold on;
    
    x_begin = poses(u).x;
    y_begin = poses(u).y;
    x_end   = poses(u+1).x;
    y_end   = poses(u+1).y;
    plot([x_begin x_end], [y_begin y_end], "r", "linewidth", 2);
    hold on;
  endfor
endif


#construct initial state
x_pose = v2t([poses(1).x; poses(1).y; poses(1).theta]);
for i=2:length(poses)
  x_pose(:,:,end+1)=v2t([poses(i).x; poses(i).y; poses(i).theta]);
endfor

#least squares loop
num_iterations = 10;
damping = 0.05;
kernel_threshold = 3.0;
[XR, XL, chi_stats, num_inliers]=do_ls(x_pose, x_lm', id_to_guess, guess_to_id,observations, transitions, length(poses), number_of_lm, num_iterations, damping, kernel_threshold);

#draw solution
figure(3)
drawLandmarks(lm_real,'g', 'fill')
lm_sol = landmark(guess_to_id(1), XL(:,1));
for u = 2:rows(XL');
  lm_sol(end+1) = landmark(guess_to_id(u), XL(:,u));
endfor
drawLandmarks(lm_guess, 'r', 'fill')
drawLandmarks(lm_sol,'b', 'fill')

#draw optimized trajectory
figure(4)
for u = 1:length(poses_real)-1
  #draw line from begin to end
  x_begin = poses_real(u).x;
  y_begin = poses_real(u).y;
  x_end   = poses_real(u+1).x;
  y_end   = poses_real(u+1).y;
  plot([x_begin x_end], [y_begin y_end], "g", "linewidth", 2);
  hold on;
  
  pose1 = t2v (XR(:,:,u));
  pose2 = t2v (XR(:,:,u+1));
  plot([pose1(1) pose2(1)], [pose1(2) pose2(2)], "r", "linewidth", 2);
  hold on;
endfor