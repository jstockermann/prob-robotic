function [initial_guess, id_to_guess, guess_to_id]  = init(poses, observations)

  #set initial pose at the origin - we don't know the map and neither our location
  #mu_t = [0;  #x coordinate
  #     0;  #y coordinate
  #    0]; #orientation theta (yaw angle)
  #mu=mu_t;  
      
  all_obs = [];

  id_to_obs_map = ones(1000, 1)*-1; #position can be found at id_to_obs_map(id+1)
  obs_to_id_map = ones(1000, 1)*-1; #obs_to_id_map(pos) gives id, pos satrs form 1
  counter=0;

  for t = 1:length(observations)

    #obtain current transition
    #transition = transitions(t);
    
    #chain up transitions
    #mu_t = transition_model(mu_t, transition.v);
    #mu(1:3,end+1) = mu_t;
    
    #obtain current observations
    observations_t = observations(t);
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
        all_obs(end+1, :) = ones(1,length(observations))*-10; #add new line to all_obs
      endif
      all_obs(id_to_obs_map(measurement.id), t) = measurement.bearing;
    endfor    
  endfor

  #parse all_obs -> triangulate between all possible measurements
  initial_guess = [];
  id_to_guess = ones(1000,1)*-1;
  guess_to_id = ones(1000,1)*-1;
  
  for i=1:rows(all_obs)
    all_tri_points=[];
    not_valid=[];
    total_measurements=0;
    for k=1:length(observations)
      if(all_obs(i,k)==-10)
        continue
      else
        total_measurements++;
        for l=k+1:length(observations)
          if(all_obs(i,l)==-10)
            continue
          else
            [lm_pos, valid] = triangulate([poses(k+1).x; poses(k+1).y], [poses(l+1).x; poses(l+1).y],...
               poses(k+1).theta + all_obs(i,k), poses(k+1).theta + all_obs(i,l));
            if(valid)
              all_tri_points(end+1,:)=lm_pos;
            else
              not_valid(end+1,:)=lm_pos;
            endif
          endif
        endfor
      endif
    endfor
    
    # mean over all possible points
    
    #printf("%i valid and %i not-valid triangulations for landmark %i; %i measurements \n" ,...
     #  rows(all_tri_points), rows(not_valid), obs_to_id_map(i), total_measurements)
    if(rows(all_tri_points)>0)
      initial_guess(end+1,:) = mean(all_tri_points,1)';
      
      guess_to_id(rows(initial_guess)) = obs_to_id_map(i);
      id_to_guess(obs_to_id_map(i)) = rows(initial_guess);
      
#    elseif(rows(not_valid)>0)
#      initial_guess(end+1) = mean(not_valid,1)';
#      
#      guess_to_id(rows(initial_guess)) = obs_to_id_map(i);
#      id_to_guess(obs_to_id_map(i)) = rows(initial_guess);
    endif  
    
    #inital_guess stays empty if only one measurement of landmark exists
    
  endfor  