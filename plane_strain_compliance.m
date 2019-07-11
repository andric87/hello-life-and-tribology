% ---- Function plane_strain_comliance computes
%------ compliance tensor for plane strain approximation
%---- As input takes compliance tensor S

function [S_plane] = plane_strain_compliance(S)

%----- Plane - strain approximation ------
for k=1:6
    for l=1:6
        S_plane(k,l)=S(k,l)-(S(k,3)*S(l,3))/S(3,3);
    end
end

% ---- Compliance tensor in plane strane approximation -----
S_plane(:,3)=[];
S_plane(3,:)=[];

