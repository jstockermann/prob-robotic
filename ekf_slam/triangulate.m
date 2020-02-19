function [lm_pos] = triangulate(mu1, mu2, theta1, theta2)


  c        = cos(theta1);
  s        = sin(theta1);
  R        = [c -s; s c];
  n1 = R*[1;0]; #unit vector pointing to direction of measurement1
  
  c        = cos(theta2);
  s        = sin(theta2);
  R        = [c -s; s c];
  n2 = R*[1;0];
  
  #solve mu1+k*n1 = mu2+l*n2
  mu = mu2-mu1;
  n = [n1,-n2];
  res = n\mu;
  
  lm_pos=res(1)*n1 + mu1;
  