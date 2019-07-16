%------Function lambda_theta rotates lambda depending
%------ of the slip plane inclination angle
function [lambda_theta] = rotate_m(lambda,theta)
c=cosd(theta);
s=sind(theta);
R = [c s 0;-1*s c 0;0 0 1];
lambda_theta=R*lambda*transpose(R);

