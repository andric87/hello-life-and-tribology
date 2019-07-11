%------Function lambda_theta rotates lambda depending
%------ of the slip plane inclination angle
function [lambda_theta] = rotate_m(lambda,theta)
cc=cosd(theta);
ss=sind(theta);
R = [cc ss 0;-1*ss cc 0;0 0 1];
lambda_theta=R*lambda*transpose(R);

