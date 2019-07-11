function [Q]=rotation_matrix(n1,n2,n3)

% --- initial basis set -----
ee_ini = [1 0 0;0 1 0; 0 0 1];

% --- rotated basis set ----
n1 = [-1 -1 2];
n2 = [1 1 1];
n3 = [-1 1 0];

ee = [n1;n2;n3];

ee_prim = [ee(1,:)/sqrt(ee(1,1)^2+ee(1,2)^2+ee(1,3)^2);...
			ee(2,:)/sqrt(ee(2,1)^2+ee(2,2)^2+ee(2,3)^2);
			ee(3,:)/sqrt(ee(3,1)^2+ee(3,2)^2+ee(3,3)^2)];

% ---- Computation of rotation matrix ---
Q = zeros(3,3);
for ii =1:3
	for jj = 1:3
		Q(ii,jj) = Q(ii,jj) + dot(ee_prim(ii,:),ee_ini(jj,:));
	end
end

end

