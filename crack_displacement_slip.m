function [u_s]=crack_displacement_slip(K_I,x,y,s1,s2,p1,p2,q1,q2,theta)

r=sqrt(x*x+y*y);
theta_atom=atan2d(y,x);

u_x = K_I*sqrt(2*r/pi)*real((1/(s1-s2))*(s1*p2*sqrt(cosd(theta_atom)+s2*sind(theta_atom))-s2*p1*sqrt(cosd(theta_atom)+s1*sind(theta_atom))) );
u_y = K_I*sqrt(2*r/pi)*real((1/(s1-s2))*(s1*q2*sqrt(cosd(theta_atom)+s2*sind(theta_atom))-s2*q1*sqrt(cosd(theta_atom)+s1*sind(theta_atom))) );

u_s = u_x*cosd(theta)+u_y*sind(theta);

end
