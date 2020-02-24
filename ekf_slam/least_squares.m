close all
clear
clc

#load dependencies
addpath "../"
addpath "../tools/g2o_wrapper"
addpath "../tools/visualization"
source "../tools/utilities/geometry_helpers_2d.m"

[_, poses, transitions, observations] = loadG2o("../datasets/BearingOnlySLAM/slam2D_bearing_only_initial_guess.g2o");

#initalisation

#set initial pose at the origin - we don't know the map and neither our location
mu_t = [0;  #x coordinate
      0;  #y coordinate
      0]; #orientation theta (yaw angle)
mu=mu_t;      
all_obs = ones(1,length(transitions))*-10; 

id_to_obs_map = ones(1000, 1)*-1;
obs_to_id_map = ones(1000, 1)*-1;  
counter=1;

for t = 1:length(transitions)

  #obtain current transition
  transition = transitions(t);
  
  #obtain current observations
  observations_t = observations(t);

  #chain up transitions
  mu_t = transition_model(mu_t, transition.v);
  mu(1:3,end+1) = mu_t;
  
  M = length(observations_t.observation);

  #if I've seen no landmarks, i do nothing
  if (M == 0)
    continue;
  endif
  
  for i=1:M
    measurement = observations_t.observation(i);
    if(id_to_obs_map(measurement.id+1)==-1) #first occurrence of id
      id_to_obs_map(measurement.id+1)=counter;
      obs_to_id_map(counter)=measurement.id+1;
      counter++;
      all_obs(end+1, :) = ones(1,length(transitions))*-10; 
      all_obs(end, i) = measurement.bearing;
    else
      all_obs(id_to_obs_map(measurement.id+1), t) = measurement.bearing;
    endif  
  endfor    
endfor

#parse all_obs
for i=4#:rows(all_obs)
  all_tri_points=[];
  for k=1:length(transitions)
    if(all_obs(i,k)==-10)
      continue
    else
      for l=k:length(transitions)
        if(all_obs(i,l)==-10)
          continue
        else
          printf("k: %i l: %i", k, l)
          [lm_pos, valid] = triangulate(mu(1:2,k+1), mu(1:2,l+1), mu(3,k+1)+all_obs(i,k), mu(3,l+1)+all_obs(i,l));
          if(valid)
            all_tri_points(end+1,:)=lm_pos;
          endif
        endif
      endfor
    endif
  endfor
  
  #mean over all possible points
  disp(all_tri_points)
  
endfor  
  