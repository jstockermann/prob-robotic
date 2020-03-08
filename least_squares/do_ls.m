# error and jacobian of a measured landmark
# input:
#   Xr: the robot pose (3x3 homogeneous matrix)
#   Xl: the landmark pose (2x1 vector, 2d pose in world frame)
#   z:  measured bearing of landmark
# output:
#   e: 1x1 is the difference between prediction and measurement
#   Jr: 1x3 derivative w.r.t a the error and a perturbation on the
#       pose
#   Jl: 1x2 derivative w.r.t a the error and a perturbation on the
#       landmark
function [e,Jr,Jl]=errorAndJacobianPoseLM(Xr,Xl,z)
  R=Xr(1:2,1:2);
  t=Xr(1:2,3);
  delta_t=Xl-t;
  meas = R*delta_t;
  z_hat = atan2(meas(2), meas(1));
  e = normalizeAngle(z_hat-z, 0);
  
  #jacobian piece w.r.t. robot
  Jr=zeros(1,3);
  faktor = 1/(delta_t(1)^2+delta_t(2)^2);
  Jr(1,1) = delta_t(2)*faktor; #de/dxr
  Jr(1,2) = -delta_t(1)*faktor; #de/dyr
  Jr(1,3) = -1; #de/dthetar

  #jacobian piece w.r.t. landmark
  Jl=zeros(1,2);
  Jl(1,1) = -delta_t(2)*faktor;
  Jl(1,2) = delta_t(1)*faktor;
endfunction;

function [e,Jr1,Jr2]=errorAndJacobianPosePose(Xr1,Xr2,trans)
  xr1 = t2v(Xr1);
  xr2 = t2v(Xr2);
  
  #xr2_pred = transition_model(xr1, trans);
  c = cos(xr1(3));
  s = sin(xr1(3));
  xr2_pred=xr2;
  ux = trans(1);
  ut = trans(3);
  xr2_pred(1) = xr1(1) + ux*c;
  xr2_pred(2) = xr1(2) + ux*s;
  xr2_pred(3) = xr1(3) + ut;
  
  f_pos = 1;
  f_ang = 1;
  delta_x = xr2_pred-xr2;
  e = f_pos*norm(delta_x(1:2)) + f_ang*abs(normalizeAngle(delta_x(3),0));
  
  Jr1 = zeros(1,3);
  Jr2 = zeros(1,3);
  
  factor = (1/(sqrt(delta_x(1)^2 + delta_x(2)^2)));
  Jr1(1) = f_pos*(delta_x(1))*factor;
  Jr1(2) = f_pos*(delta_x(2))*factor;
  Jr1(3) = 2*f_pos*(ux*c*delta_x(2)-ux*s*delta_x(1))*factor + f_ang*sign(delta_x(3));
  
  Jr2(1) = -Jr1(1);
  Jr2(2) = -Jr1(2);
  Jr2(3) = -f_ang*sign(delta_x(3));

endfunction

function [e,Jr1,Jr2]=errorAndJacobianPosePose_flatten(Xr1,Xr2,trans)
  xr1 = t2v(Xr1);
  xr2 = t2v(Xr2);
  
  #xr2_pred = transition_model(xr1, trans);
  c = cos(xr1(3));
  s = sin(xr1(3));
  xr2_pred=xr2;
  ux = trans(1);
  ut = trans(3);
  xr2_pred(1) = xr1(1) + ux*c;
  xr2_pred(2) = xr1(2) + ux*s;
  xr2_pred(3) = xr1(3) + ut;
  Xr2_pred = v2t(xr2_pred);
  
  e = reshape(Xr2_pred(1:2,1:3),6,1) - reshape(Xr2(1:2,1:3),6,1);
  
  Jr1 = zeros(6,3);
  Jr2 = zeros(6,3);
  
  Jr1(1,3) = -sin(xr1(3) + ut);
  Jr1(2,3) = cos(xr1(3) + ut);
  Jr1(3,3) = -cos(xr1(3) + ut);
  Jr1(4,3) = -sin(xr1(3) + ut);
  Jr1(5,:) = [1, 0, -ux * sin(xr1(1))];
  Jr1(6,:) = [0, 1, ux * cos(xr1(1))];
  
  Jr2(1,3) = sin(xr2(3));
  Jr2(2,3) = -cos(xr2(3));
  Jr2(3,3) = cos(xr2(3));
  Jr2(4,3) = sin(xr2(3));
  Jr2(5,:) = [-1, 0, 0];
  Jr2(6,:) = [0, -1, 0];
  

endfunction


# implementation of the boxplus
# applies a perturbation to a set of landmarks and robot poses
# input:
#   XR: the robot poses (4x4xnum_poses: array of homogeneous matrices)
#   XL: the landmark pose (3xnum_landmarks matrix of landmarks)
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   dx: the perturbation vector of appropriate dimensions
#       the poses come first, then the landmarks
# output:
#   XR: the robot poses obtained by applying the perturbation
#   XL: the landmarks obtained by applying the perturbation
function [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx)
  pose_dim = 3;
  landmark_dim = 2;
  for(pose_index=1:num_poses)
    pose_matrix_index=1+(pose_index-1)*pose_dim;
    dxr=dx(pose_matrix_index:pose_matrix_index+pose_dim-1);
    XR(:,:,pose_index)=v2t(0.2*dxr)*XR(:,:,pose_index);
  endfor;
  for(landmark_index=1:num_landmarks)
    landmark_matrix_index=1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;
    dxl=dx(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,:);
    XL(:,landmark_index)+=dxl;
  endfor;
endfunction;



# implementation of the optimization loop with robust kernel
# applies a perturbation to a set of landmarks and robot poses
# input:
#   XR: the initial robot poses (3x3xnum_poses: array of homogeneous matrices)
#   XL: the initial landmark estimates (2xnum_landmarks matrix of landmarks)
#   Z:  the measurements (3xnum_measurements)
#   associations: 2xnum_measurements. 
#                 associations(:,k)=[p_idx,l_idx]' means the kth measurement
#                 refers to an observation made from pose p_idx, that
#                 observed landmark l_idx
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   num_iterations: the number of iterations of least squares
#   damping:      damping factor (in case system not spd)
#   kernel_threshod: robust kernel threshold

# output:
#   XR: the robot poses after optimization
#   XL: the landmarks after optimization
#   chi_stats: array 1:num_iterations, containing evolution of chi2
#   num_inliers: array 1:num_iterations, containing evolution of inliers
function [XR, XL, chi_stats, num_inliers]=do_ls(XR, XL, id_to_guess, guess_to_id, observations, transitions,
							num_poses, 
							num_landmarks, 
							num_iterations, 
							damping, 
							kernel_threshold)
  pose_dim = 3;
  landmark_dim = 2;

  chi_stats=zeros(1,num_iterations);
  num_inliers=zeros(1,num_iterations);
  # size of the linear system
  system_size=pose_dim*num_poses+landmark_dim*num_landmarks; 
  for (iteration=1:num_iterations)
    H=zeros(system_size, system_size);
    b=zeros(system_size,1);
    chi_stats(iteration)=0;
    for (t = 1:length(observations))
      observation_t=observations(t);
      for (obs = 1: length(observation_t))
        observation = observation_t.observation(obs);
        pose_index = t+1;
        if(observation.id == 0 || id_to_guess(observation.id)==-1) 
          continue;
        endif
        #val = debug_on_error()
        landmark_index = id_to_guess(observation.id);
        z=observation.bearing;
        Xr=XR(:,:,pose_index);
        Xl=XL(:,landmark_index);
        [e,Jr,Jl] = errorAndJacobianPoseLM(Xr, Xl, z);
        chi=e'*e;
        if (chi>kernel_threshold)
          e*=sqrt(kernel_threshold/chi);
          chi=kernel_threshold;
        else
          num_inliers(iteration)++;
        endif;
        chi_stats(iteration)+=chi;

        Hrr=Jr'*Jr;
        Hrl=Jr'*Jl;
        Hll=Jl'*Jl;
        br=Jr'*e;
        bl=Jl'*e;

        pose_matrix_index = 1 + (pose_index-1)*pose_dim;
        landmark_matrix_index = 1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;

        H(pose_matrix_index:pose_matrix_index+pose_dim-1,
    pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hrr;

        H(pose_matrix_index:pose_matrix_index+pose_dim-1,
    landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hrl;

        H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
    landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hll;

        H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
    pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hrl';

        b(pose_matrix_index:pose_matrix_index+pose_dim-1)+=br;
        b(landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=bl;
      endfor
    endfor
    for(p = 1:(num_poses-1))
      Xr1=XR(:,:,p);
      Xr2=XR(:,:,p+1);
      trans = transitions(p).v; 
      #[e,Jr1,Jr2]=errorAndJacobianPosePose(Xr1,Xr2,trans);
      [e,Jr1,Jr2]=errorAndJacobianPosePose_flatten(Xr1,Xr2,trans);
      chi=e'*e;
      if (chi>kernel_threshold)
        e*=sqrt(kernel_threshold/chi);
        chi=kernel_threshold;
      else
        num_inliers(iteration)++;
      endif;
      chi_stats(iteration)+=chi;
      
      Hrr=Jr1'*Jr1;
      Hrl=Jr1'*Jr2;
      Hll=Jr2'*Jr2;
      br=Jr1'*e;
      bl=Jr2'*e;

      pose1_matrix_index = 1 + (p-1)*pose_dim;
      pose2_matrix_index = 1 + (p)*pose_dim;
      landmark_matrix_index = 1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;

      H(pose1_matrix_index:pose1_matrix_index+pose_dim-1,
  pose1_matrix_index:pose1_matrix_index+pose_dim-1)+=Hrr;

      H(pose1_matrix_index:pose1_matrix_index+pose_dim-1,
  pose2_matrix_index:pose2_matrix_index+pose_dim-1)+=Hrl;

      H(pose2_matrix_index:pose2_matrix_index+pose_dim-1,
  pose2_matrix_index:pose2_matrix_index+pose_dim-1)+=Hll;

      H(pose2_matrix_index:pose2_matrix_index+pose_dim-1,
  pose1_matrix_index:pose1_matrix_index+pose_dim-1)+=Hrl';

      b(pose1_matrix_index:pose1_matrix_index+pose_dim-1)+=br;
      b(pose2_matrix_index:pose2_matrix_index+pose_dim-1)+=bl;
    endfor
    
    H+=eye(system_size)*damping;
    dx=zeros(system_size,1);

    % we solve the linear system, blocking the first pose
    % this corresponds to "remove" from H and b the locks
    % of the 1st pose, while solving the system

    dx(pose_dim+1:end)=-(H(pose_dim+1:end,pose_dim+1:end)\b(pose_dim+1:end,1));
    [XR, XL]=boxPlus(XR,XL,num_poses, num_landmarks, dx);
    
    figure(4)
    clf(4)
    hold off;
    for u = 1:num_poses-1
      #draw line from begin to end     
      
      pose1 = t2v (XR(:,:,u));
      pose2 = t2v (XR(:,:,u+1));
      plot([pose1(1) pose2(1)], [pose1(2) pose2(2)], "r", "linewidth", 2);
      hold on;
    endfor
    drawnow();
    pause(1);
    chi_stats(iteration)
  endfor
endfunction