%---Function stroh_tensor computes Stroh energy tensor
%-- based on fundamental elasticity tensor. For more theoretical detals see 
%--- "Anisotropic elasticity" Ting
%--- As input takes material stiffness tensor C (6x6) in the coordinates of interest 
% ---------------------------------------------------------------
function [lambda]=stroh_tensor(C)

% ---- Finding fundamental elasticity matrix N------
Q=[C(1,1) C(1,6) C(1,5);C(1,6) C(6,6) C(5,6);C(1,5) C(5,6) C(5,5)];
R=[C(1,6) C(1,2) C(1,4);C(6,6) C(2,6) C(4,6);C(5,6) C(2,5) C(4,5)];
T=[C(6,6) C(2,6) C(4,6);C(2,6) C(2,2) C(2,4);C(4,6) C(2,4) C(4,4)];
N1=-1*T^(-1)*transpose(R);
N2=T^(-1);
N3=R*T^(-1)*transpose(R)-Q;
N=[N1 N2;N3 transpose(N1)];

[u,v] = eig(N); % finding eigen values of the fundamental elasticity matrix

A = [];
B = [];
p = [];

for jj = 1:2:5
    a = [];
    b = [];
    for ii = 1:3
	a1 = u(ii,jj);
	b1 = u(ii+3,jj);
	a = [a;a1];
	b = [b;b1];
    end

    A = [A,a];
    B = [B,b];

end
A
B
L = real(i*A*B^(-1));
lambda = 0.5*L;

end


