close all
clear
clc

#load dependencies
addpath "../"
addpath "../tools/g2o_wrapper"
addpath "../tools/visualization"
addpath "../ekf_slam"
source "../tools/utilities/geometry_helpers_2d.m"

[_, poses, transitions, observations] = loadG2o("../datasets/BearingOnlySLAM/slam2D_bearing_only_initial_guess.g2o");

#initalisation

inital_guess = init(poses, observations);
  