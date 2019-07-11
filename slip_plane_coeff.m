function [F12]=slip_plane_coeff(s1, s2, theta)

c=cosd(theta);
s=sind(theta);

F1=real(((s1*s2)/(s1-s2))*(s2/sqrt(c+s2*s)-s1/sqrt(c+s1*s)));
F2=real((1/(s1-s2))*(s1/sqrt(c+s2*s)-s2/sqrt(c+s1*s)));
F3=real(((s1*s2)/(s1-s2))*(1/sqrt(c+s1*s)-1/sqrt(c+s2*s)));
F12=c*(c*F3-s*F1)+s*(c*F2-s*F3);
end
