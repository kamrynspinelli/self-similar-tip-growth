function [X_dot] = solver_parabolic_degenerate_fast(X,LUT,L0,R0,K,MU,ext_verts,ext_force_status)
global P 
%X is a row vector
%X_dot is a column vector
N = size(LUT,2); % Number of vertices
rv = reshape(X',N,2)'; %Vertex positions in the usual format 2 by N
arcs = ParabolicArc.all_degenerate_arcs(rv, K, MU, L0, R0);
r_dot = zeros(2,N); % preallocate matrix of force vectors

for i=1:N % for each marker point,
    % add up all the forces on that point
    r_dot(:,i) = ParabolicArc.Tvec2D(rv(:,i), arcs(i), arcs(i+1)) ...
        + ParabolicArc.Tvec3DrR(rv(:,i), arcs(i), arcs(i+1)) ...
        + ParabolicArc.Tvec3DnR(rv(:,i), arcs(i), arcs(i+1)) ...
        + ParabolicArc.Pvec(rv(:,i), arcs(i), arcs(i+1), P);
end

r_dot(2,end) = 0;
r_dot(1,end) = 10*(P/2*arcs(end-1).b(2)-arcs(end-1).sigmaS*cos(arcs(end-1).beta-pi/2)); % define tip behavior to match analytic linear segments model
if ext_force_status == 1
    r_dot = r_dot + external_force(rv);
end
% disp(norm(r_dot, Inf));
% if norm(r_dot, Inf) > 20
%     disp("big")
% end

X_dot = reshape(r_dot',2*N,1);

end