% ---Function crack_anisotropic_input computes
%----input for crack stress and diplacement field

function [s1,s2,p1,p2,q1,q2]=crack_anisotropic_input(S_plane)
format long
%syms x
%Eq = S_plane(1,1)*x^4-2*S_plane(1,5)*x^3+(2*S_plane(1,2)+S_plane(5,5))*x^2-2*S_plane(2,5)*x+S_plane(2,2);
%sol = double(solve(Eq,x));
sol = roots([S_plane(1,1) -2*S_plane(1,5) 2*S_plane(1,2)+S_plane(5,5) -2*S_plane(2,5) S_plane(2,2)]);
   if imag(sol(1,1))>0
	 s1=sol(1,1);
   else s1=sol(2,1);
   end
   if imag(sol(3,1))>0
	 s2=sol(3,1);
   else s2=sol(4,1);
   end

% ----- Computing parameters for anisotropic displacement field around the crack
p1 = S_plane(1,1)*s1^2+S_plane(1,2)-S_plane(1,5)*s1;
p2 = S_plane(1,1)*s2^2+S_plane(1,2)-S_plane(1,5)*s2;
q1 = S_plane(1,2)*s1+S_plane(2,2)/s1-S_plane(2,5);
q2 = S_plane(1,2)*s2+S_plane(2,2)/s2-S_plane(2,5);
format
end
