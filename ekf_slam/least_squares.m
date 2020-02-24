close all
#clear
#clc

#load dependencies
addpath "../"
addpath "../tools/g2o_wrapper"
addpath "../tools/visualization"
source "../tools/utilities/geometry_helpers_2d.m"

#[_, poses, transitions, observations] = loadG2o("../datasets/BearingOnlySLAM/slam2D_bearing_only_initial_guess.g2o");

#initalisation

#set initial pose at the origin - we don't know the map and neither our location
#mu_t = [0;  #x coordinate
 #     0;  #y coordinate
  #    0]; #orientation theta (yaw angle)
#mu=mu_t;  
mu=[poses(1).x; poses(1).y; poses(1).theta];
    
all_obs = []; 
all_obs_c = cell(1);

id_to_obs_map = ones(1000, 1)*-1; #position can be found at id_to_obs_map(id+1)
obs_to_id_map = ones(1000, 1)*-1; #obs_to_id_map(pos) gives id, pos satrs form 1
counter=0;

for t = 1:length(transitions)

  #obtain current transition
  transition = transitions(t);
  
  #obtain current observations
  observations_t = observations(t);

  #chain up transitions
  #mu_t = transition_model(mu_t, transition.v);
  #mu(1:3,end+1) = mu_t;
  
  #use poses
  mu(1:3,end+1) = [poses(t+1).x; poses(t+1).y; poses(t+1).theta];
  
  M = length(observations_t.observation);

  #if I've seen no landmarks, i do nothing
  if (M == 0)
    continue;
  endif
  
  for i=1:M
    measurement = observations_t.observation(i);
    if (measurement.id <= 0) #ignore lm with zero id
      continue;
    endif
    if(id_to_obs_map(measurement.id)==-1) #first occurrence of id
      counter++;
      id_to_obs_map(measurement.id)=counter;
      obs_to_id_map(counter)=measurement.id;      
      all_obs(end+1, :) = ones(1,length(transitions))*-10; #add new line to all_obs
    endif
    all_obs(id_to_obs_map(measurement.id), t) = measurement.bearing;
  endfor    
endfor

#parse all_obs
inital_guess = cell(1);
for i=1:rows(all_obs)
  all_tri_points=[];
  not_valid=[];
  total_measurements=0;
  for k=1:length(transitions)
    if(all_obs(i,k)==-10)
      continue
    else
      total_measurements++;
      for l=k+1:length(transitions)
        if(all_obs(i,l)==-10)
          continue
        else
          # printf("k: %i l: %i /n" , k, l)
          [lm_pos, valid] = triangulate(mu(1:2,k+1), mu(1:2,l+1), mu(3,k+1)+all_obs(i,k), mu(3,l+1)+all_obs(i,l));
          if(valid)
            all_tri_points(end+1,:)=lm_pos;
          else
            not_valid(end+1,:)=lm_pos;
          endif
        endif
      endfor
    endif
  endfor
  
  #mean over all possible points
  printf("%i valid and %i not-valid triangulations for landmark %i; %i measurements \n" ,...
     rows(all_tri_points), rows(not_valid), obs_to_id_map(i), total_measurements)
  if(rows(all_tri_points)>0)
    inital_guess(obs_to_id_map(i)) = mean(all_tri_points,1);
  elseif(rows(not_valid)>0)
    inital_guess(obs_to_id_map(i)) = mean(not_valid,1);
  endif  
  
  #inital_guess stays empty if only one measurement of landmark exists
  
endfor  
  