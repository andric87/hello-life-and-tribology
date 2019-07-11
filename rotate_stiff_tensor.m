function [C_r]=rotate_stiff_tensor(C,m1,m2,m3)

% --- Crystal orinentation in cubic coordinates
e1=transpose([1 0 0]);
e2=transpose([0 1 0]);
e3=transpose([0 0 1]);

%----Desired crystal orientation ------
m11=m1/sqrt(m1(1,1)^2 + m1(1,2)^2 + m1(1,3)^2); %lattice orientation
m22=m2/sqrt(m2(1,1)^2 + m2(1,2)^2 + m2(1,3)^2); %lattice orientation
m33=m3/sqrt(m3(1,1)^2 + m3(1,2)^2 + m3(1,3)^2); %lattice orientaion

%-----Rotation matrix
Q=[m11*e1 m11*e2 m11*e3;m22*e1 m22*e2 m22*e3;m33*e1 m33*e2 m33*e3];

K1=[Q(1,1)^2 Q(1,2)^2 Q(1,3)^2;...
	Q(2,1)^2 Q(2,2)^2 Q(2,3)^2;...
	Q(3,1)^2 Q(3,2)^2 Q(3,3)^2];
K2=[Q(1,2)*Q(1,3) Q(1,3)*Q(1,1) Q(1,1)*Q(1,2);...
	Q(2,2)*Q(2,3) Q(2,3)*Q(2,1) Q(2,1)*Q(2,2);...
	Q(3,2)*Q(3,3) Q(3,3)*Q(3,1) Q(3,1)*Q(3,2)];
K3=[Q(2,1)*Q(3,1) Q(2,2)*Q(3,2) Q(2,3)*Q(3,3);...
	Q(3,1)*Q(1,1) Q(3,2)*Q(1,2) Q(3,3)*Q(1,3);...
	Q(1,1)*Q(2,1) Q(1,2)*Q(2,2) Q(1,3)*Q(2,3)];
K4=[Q(2,2)*Q(3,3)+Q(2,3)*Q(3,2) Q(2,3)*Q(3,1)+Q(2,1)*Q(3,3) Q(2,1)*Q(3,2)+Q(2,2)*Q(3,1);...
	Q(3,2)*Q(1,3)+Q(3,3)*Q(1,2) Q(3,3)*Q(1,1)+Q(3,1)*Q(1,3) Q(3,1)*Q(1,2)+Q(3,2)*Q(1,1);...
	Q(1,2)*Q(2,3)+Q(1,3)*Q(2,2) Q(1,3)*Q(2,1)+Q(1,1)*Q(2,3) Q(1,1)*Q(2,2)+Q(1,2)*Q(2,1)];

KK=[K1 2*K2;K3 K4]; %---rotation matrix for 4th order tensor rotation (given in a matrix form 6x6)
%------------------------------------------------------------------------------------------------
C_r = KK*C*transpose(KK);
end

