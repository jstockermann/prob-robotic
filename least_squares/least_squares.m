close all
#clear
#clc

#load dependencies
addpath "../"
addpath "../tools/g2o_wrapper"
addpath "../tools/visualization"
addpath "../ekf_slam"
source "../tools/utilities/geometry_helpers_2d.m"

[_, poses, transitions, observations] = loadG2o("../datasets/BearingOnlySLAM/slam2D_bearing_only_initial_guess.g2o");

#initalization

[initial_guess, id_to_guess, guess_to_id] = init(poses, observations);

#----------------------
# intialization works but least squares is just a draft
#----------------------

#construct initial state
x_lm = initial_guess;
x_pose = [];
for i=2:length(poses)
  x_pose[end+1]=[poses(i).x; poses(i).y; poses(i).theta];
endfor

#least squares loop
chi = 2;
while(chi>1)
  #pose-LM
  chi_lm=0; 
  
  dim = length(x_lm)+length(x_pose);
  H=zeros(dim,dim);  b=zeros(dim,1); %accumulators for H and b
  for(i=1:length(observations))
     observations_t = observations(i)
     [e,J]=errorAndJacobian_LM(x_lm, x_pose, observations_t); %compute e and J 
     H+=J'*J;            %assemble H and B
     b+=J'*e;
     chi_lm+=e'*e;          %update cumulative error
  endfor
  dx=­H\b;               %solve the linear system
  x_new=x+dx;            %apply update
  