function mg_fracture
clear;clc;
format long

% ---- Material properties ----
% ----- Mg sstiffness tensor entry (Transveresly isotropic)
% ----- Units: (GPa) !!!!!!!!!    ----
C11 = 68.5784405816146; 
C22 = C11;
C12 = 29.6300165274605; 
C13 = 22.5287657640503;
C23 = C13;
C33 = 65.5535542657224;
C44 = 19.1515545182744;
C55 = C44;
C66 = (C11-C12)/2;

% ---- MG stiffness tensor in Voigt notation
C = [C11 C12 C13 0 0 0;...
     C12 C22 C23 0 0 0;...
     C13 C23 C33 0 0 0;...
     0 0 0 C44 0 0;...
     0 0 0 0 C44 0;...
     0 0 0 0 0 C66];

% ---- Basis set depedning of the crystal orientation ----
n1 = [1 0 0]; % --- what I need in the simulation
n2 = [0 0 1];
n3 = [0 1 0];
theta_slip = 61.916;
phi_1 = 15.205;
phi_2 = phi_1+60;

% --- surface and stacking fault energies ----
gamma_s = 0.611;

gamma_usf = 0.2042
gamma_se = 0.449

% --- computin rotation matrix ----
C_rot = rotate_stiff_tensor(C,n1,n2,n3);

S_rot = inv(C_rot);

S_plane = plane_strain_compliance(S_rot);

[s1 s2 p1 p2 q1 q2] = crack_anisotropic_input(S_plane); % --- parameter for crack-tip stress and displacement fields
% --- printing input parameters for fracture simulations 
fileID = fopen('fracture_input_Mg.txt','w');
fprintf(fileID,'s1=%.6f %.6f*i\n',real(s1),imag(s1));
fprintf(fileID,'s2=%.6f %.6f*i \n',real(s2),imag(s2));
fprintf(fileID,'p1=%.6f %.6f*i \n',real(p1),imag(p1));
fprintf(fileID,'p2=%.6f %.6f*i \n',real(p2),imag(p2));
fprintf(fileID,'q1=%.6f %.6f*i\n',real(q1),imag(q1));
fprintf(fileID,'q2=%.6f %.6f*i\n',real(q2),imag(q2));
fclose(fileID);

F12 = slip_plane_coeff(s1, s2, theta_slip); % ---- correction for slip plane

lambda = stroh_tensor_B(C_rot); %---Stroh energy tensor
inv_lambda = inv(lambda);
lambda_theta = rotate_m(lambda,theta_slip);

s_phi1 = [cosd(phi_1) 0 sind(phi_1)]; % DIrection of the Burgers vector for the 1st emission
s_phi2 = [cosd(phi_2) 0 sind(phi_2)]; % DIrection of the Burgers vector for the 2nd emission

o_first = s_phi1*lambda_theta^(-1)*transpose(s_phi1);
o_second = s_phi2*lambda_theta^(-1)*transpose(s_phi2);

K_Ic = sqrt(2*(gamma_s/lambda(2,2)*10^(9)))*10^(-6) % Griffith value (MPa sqrt(m))

if gamma_se/gamma_usf>3.45
    G_Ie = (0.145*gamma_se+0.5*gamma_usf);
else
    G_Ie = gamma_usf;
end

% Sttress intensity for  emission
K_Ie = (sqrt(G_Ie*o_first*10^9)*10^(-6))/F12/cosd(phi_1)
K_Ie_Rice = (sqrt(gamma_usf*o_first*10^9)*10^(-6))/F12/cosd(phi_1)

gamma_se/gamma_usf

