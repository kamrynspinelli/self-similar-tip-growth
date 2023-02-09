function [X_dot] = solver_fast(X,LUT,L0,R0,K,MU,ext_verts,ext_force_status)
global P 
angle_0=asin(2*R0(end)/L0(end));
M = size(LUT,1);
N = size(LUT,2); % Number of vertices
ry = ones(M,2);
ry = [ry(:,1).*0 ry(:,2)];
rv = reshape(X',N,2)'; %Vertex positions in the usual format 2 by N
rb = LUT*rv';%boundary vectors
LUTP = LUT>0;
LUTM = LUT<0;
rb1 = (LUTP+LUTM)*rv(2,:)';%yL+yR 
rbL = LUTP*rv(2,:)';
rbR = LUTM*rv(2,:)';
Nb = size(ext_verts,2);% Number of boundary vertices
% Generate boundary vertice tangential movement matrix
Mb = zeros(2,Nb);
Mb(1,rv(1,ext_verts)>1e-10) = 1;
Mb(2,rv(2,ext_verts)>1e-10) = 1;
%Compute instantaneous tension vectors
D = sqrt(sum(rb.^2,2)); %Vector of edge lengths.
rb(:,1) = rb(:,1)./D;
rb(:,2) = rb(:,2)./D;%director
rbn = [rb(:,2) -rb(:,1)];%outward normal
rm = 0.5*rb1;
angle = -atan(rb(:,2)./rb(:,1));
lnL = log(rbR./rm);
lnR = log(rm./rbL);
%Fl = K.*((D./L0).^2-(D./L0.*rm./R0).^(-4));
%Fr = K.*((rm./R0).^2-(D./L0.*rm./R0).^(-4));
Fl = K.*((D./L0.*rm./R0)-1)+0.5*MU.*((R0.^2./rm.^2)-L0.^2./D.^2);
Fr = K.*((D./L0.*rm./R0)-1)+0.5*MU.*(L0.^2./D.^2-(R0.^2./rm.^2));

%The meridional tension in 2D
rb2D = [Fl Fl].*(-rb);
%rb2D = (-rb);
Tvec2D = rb2D';

%The 3D effect of tension anisotropy
%rb3DrL = [1/2*D./rm.*(Fl-Fr) 1/2*D./rm.*(Fl-Fr)].*(ry);
% A=1/2*D./rm
% B=-lnL./sin(angle) 
% C=-lnR./sin(angle) 
rb3DrL = [lnL.*(Fl-Fr)./sin(-angle) lnL.*(Fl-Fr)./sin(-angle)].*(ry);
%rb3DtL = [lnL lnL].*(rb);
rb3DrL(end,:)=0;%remove the singularity at r=0
Tvec3DrL = rb3DrL';
%rb3DrR = [1/2*D./rm.*(Fl-Fr) 1/2*D./rm.*(Fl-Fr)].*(ry);
rb3DrR = [lnR.*(Fl-Fr)./sin(-angle) lnR.*(Fl-Fr)./sin(-angle)].*(ry);
%rb3DtR = [lnR lnR].*(rb);
Tvec3DrR = rb3DrR';

%The 3D effect of meridional tension
%rb3DnL = [1/2*D./rm.*Fl.*cos(angle) 1/2*D./rm.*Fl.*cos(angle)].*(-rbn);
rb3DnL = [lnL.*Fl.*cos(-angle)./sin(-angle) lnL.*Fl.*cos(-angle)./sin(-angle)].*(-rbn);
%rb3DrL = [lnL./sin(angle) lnL./sin(angle)].*(-ry);
rb3DnL(end,:)=0;%remove the singularity at r=0
Tvec3DnL = rb3DnL';

%rb3DnR = [1/2*D./rm.*Fl.*cos(angle) 1/2*D./rm.*Fl.*cos(angle)].*(-rbn);
rb3DnR = [lnR.*Fl.*cos(-angle)./sin(-angle) lnR.*Fl.*cos(-angle)./sin(-angle)].*(-rbn);
%rb3DrR = [lnR./sin(angle) lnR./sin(angle)].*(-ry);
Tvec3DnR = rb3DnR';

%The pressure in the normal direction
Pvec = 1/2*P*[rbn(:,1).*D rbn(:,2).*D];
Pvec = Pvec';
% a = Tvec2D*LUT
% b = Pvec*abs(LUT)
% c = Tvec3DrR*LUTP
% d = Tvec3DtR*LUTP
% e = (Tvec3DrL)*LUTM
% f = Tvec3DtL*LUTM
%
r_dot = Pvec*abs(LUT)+Tvec2D*LUT+(Tvec3DnR+Tvec3DrR)*LUTP+(Tvec3DnL+Tvec3DrL)*LUTM;
%r_dot = -(Tvec2D)*LUT+(Tvec3Dn)*abs(LUT)+(Tvec3Dt)*abs(LUT)-Pvec*abs(LUT);
r_dot(:,ext_verts) = 2.*Mb.*r_dot(:,ext_verts);%Making the sliding boundary
r_dot(1,ext_verts(2)) = 10*(P/2*rm(end)-Fl(end)*cos(-angle(end)));% The tip is a cone, the boundary condition fixed, multiply by 10 to make it faster
%r_dot(1,ext_verts(2)) = (rv(2,ext_verts(2)-1)-rv(2,ext_verts(2)))-tan(-angle_0)*(rv(1,ext_verts(2)-1)-rv(1,ext_verts(2)));%The tip condition as a geometric constraint
if ext_force_status == 1
    r_dot = r_dot + external_force(rv);
end

X_dot = reshape(r_dot',2*N,1);

end