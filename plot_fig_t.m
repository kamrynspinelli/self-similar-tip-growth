function [Tmax Tmin] = plot_fig_t(rv_t,LUT_t,Tl_t,Tr_t,R0,K,mode)
%Plot only the lattice structures
global gr gl k kr kl L Fext variance_factor
%filename = ['gr_',num2str(gr),'gl_',num2str(gl),'kr_',num2str(kr/L),'kl_',num2str(kl/L)];
filename = ['k_',num2str(k),'kr_',num2str(kr/L),'kl_',num2str(kl/L),'_Fext',num2str(Fext),'_var',num2str(variance_factor),'_',mode];

%% Find image length scale (largest length scale of system).
MIN_X = Inf;
MIN_Y = Inf;
MAX_X = -Inf;
MAX_Y = -Inf;

MIN_X0 = Inf*ones(1,length(rv_t));
MIN_Y0 = Inf*ones(1,length(rv_t));
MAX_X0 = -Inf*ones(1,length(rv_t));
MAX_Y0 = -Inf*ones(1,length(rv_t));

for jj=1:length(rv_t)
    rvt = rv_t(jj).dat;
    MIN_X0(jj) = min(min(rvt(1,:)),MIN_X0(jj));
    MAX_X0(jj) = max(max(rvt(1,:)),MAX_X0(jj));
    MIN_Y0(jj) = min(min(rvt(2,:)),MIN_Y0(jj));
    MAX_Y0(jj) = max(max(rvt(2,:)),MAX_Y0(jj));
end


for jj=2:length(rv_t)
    rvt = rv_t(jj).dat;
    MIN_X = min(min(rvt(1,:)),MIN_X);
    MAX_X = max(max(rvt(1,:)),MAX_X);
    MIN_Y = min(min(rvt(2,:)),MIN_Y);
    MAX_Y = max(max(rvt(2,:)),MAX_Y);
end

LARGERDIMENSION = max(2*(MAX_Y-MIN_Y),(MAX_X-MIN_X));
MAX_X = LARGERDIMENSION+MIN_X;
MAX_Y = 1/2*LARGERDIMENSION;
MIN_Y = -1/2*LARGERDIMENSION;


Tmin = Inf;
Tmax = -Inf;

Kmin = Inf;
Kmax = -Inf;
Kmin =min(Kmin,K);
Kmax =max(Kmax,K);
KK = (K-Kmin)/(Kmax-Kmin);
nK = length(K);

for jj=1:length(LUT_t)
    Tlt = Tl_t(jj).dat;
    Trt = Tr_t(jj).dat;
    Tmin = min(min(Tlt),Tmin);
    Tmax = max(max(Tlt),Tmax);
    Tmin = min(min(Trt),Tmin);
    Tmax = max(max(Trt),Tmax);
end
   rvnew = rv_t(1).dat;
for kk = 1:size(rvnew,2)-1
   rvnew(2,kk) = (R0(kk));
end
   % rv_t(1).dat = rvnew;
    Tl_t(1).dat = Tmin*ones(length(Tl_t(1).dat),1);
    Tr_t(1).dat = Tmin*ones(length(Tr_t(1).dat),1);

v = VideoWriter(['Video_',filename,'.avi']);
open(v)

%% Plot image
for jj=1:length(rv_t)
    Kgradient = (MAX_X0(jj)-MIN_X0(jj));
    figure;
    %xlim([MIN_X MAX_X]);
    %ylim([MIN_Y MAX_Y]);
    axis equal
    axis([MIN_X MAX_X*(1+0.05) MIN_Y MAX_Y]);
    if mode == "tension"
        ylabel({'Circumferential (lower half) and meridional (upper half) tension'},'FontSize',12,'FontWeight','bold');
    elseif mode == "curvature"
        ylabel({'Mean (lower half) and Gaussian (upper half) curvature'},'FontSize',12,'FontWeight','bold');
    end
    hold on
    rvt = rv_t(jj).dat;
    LUTt = LUT_t(jj).dat;
    Tlt = Tl_t(jj).dat;
    Trt = Tr_t(jj).dat;
for ii=1:size(LUTt,1)
    index = find(LUTt(ii,:));
    x1 = rvt(1,index(1));
    y1 = rvt(2,index(1));
    x2 = rvt(1,index(2));
    y2 = rvt(2,index(2));
    h(1) = plot([x1 x2],[y1 y2]);
    h(1).LineWidth = 5;
    T = (Tlt(ii)-Tmin)/(Tmax-Tmin);
    x=ones(1,10);
    for n=1:10
        x(n)=255;
    end
    y=ones(1,10);
    for n=1:10
        y(n)=222;
    end
    z=ones(1,10);
    z(1)=230.14;
    for n=2:10
        z(n)=z(n-1)-25.571;
    end
    a=ones(1,50);
    for n=1:50
        a(n)=255;
    end
    b=ones(1,50);
    b(1)=222;
    for n=2:50
        b(n)=b(n-1)-4.5306;
    end
    c=ones(1,50);
    for n=1:50
        c(n)=0;
    end
    d=ones(1,10);
    e=ones(1,10);
    f=ones(1,10);
    d(1)=255;
    for n=2:10
        d(n)=d(n-1)-11.6;
    end
    for n=1:10
        e(n)=0;
        f(n)=0;
    end
    if T>2/3
        r = round(3*(T-2/3)/(1/9))+1;
        c1 = d(r)/255;
        c2 = e(r)/255;
        c3 = f(r)/255;
    elseif (T>1/3)&&(T<2/3)
        r = round(3*(T-1/3)/(1/49))+1;
        c1 = a(r)/255;
        c2 = b(r)/255;
        c3 = c(r)/255;
    else
        r = round(3*T/(1/9))+1;
        c1 = x(r)/255;
        c2 = y(r)/255;
        c3 = z(r)/255;
    end
    h(1).Color = [c1 c2 c3];
    %if T>0.5
    %h(1).Color = [1-T 0 0];
    %else
    %h(1).Color = [1 1-T 1-T];
    %end
   
    h(2) = plot([x1 x2],[-y1 -y2]);
    h(2).LineWidth = 5;
    T = (Trt(ii)-Tmin)/(Tmax-Tmin);
    if T>2/3
        r = round(3*(T-2/3)/(1/9))+1;
        c1 = d(r)/255;
        c2 = e(r)/255;
        c3 = f(r)/255;
    elseif (T>1/3)&&(T<2/3)
        r = round(3*(T-1/3)/(1/49))+1;
        c1 = a(r)/255;
        c2 = b(r)/255;
        c3 = c(r)/255;
    else
        r = round(3*T/(1/9))+1;
        c1 = x(r)/255;
        c2 = y(r)/255;
        c3 = z(r)/255;
    end
    h(2).Color = [c1 c2 c3];
    %h(3) = plot([0 1],[0 0]);
    %h(3).LineWidth = 5;
    %h(3).Color = [1 1 0];

    %if T>0.5
    %h(2).Color = [1-T 0 0];
    %else
    %h(2).Color = [1 1-T 1-T];
    %end
end

    for n=1:20
        h(n+2)=plot([MAX_X*1.05 MAX_X*1.05],[-0.5+(0+n-1)/60 -0.5+(n/60)]);
        h(n+2).LineWidth=10;
        c1=255;
        c2=222;
        c3=230.14;
        h(n+2).Color=[c1/255 c2/255 (c3-(n-1)*11.507)/255];
    end
    for n=1:20
        h(n+22)=plot([MAX_X*1.05 MAX_X*1.05],[-0.5+1/3+(n-1)/60 -0.5+1/3+n/60]);
        h(n+22).LineWidth=10;
        c1=255;
        c2=222;
        c3=0;
        h(n+22).Color=[c1/255 (c2-(n-1)*11.1)/255 c3/255];
    end
    for n=1:20
        h(n+42)=plot([MAX_X*1.05 MAX_X*1.05],[-0.5+2/3+(n-1)/60 -0.5+2/3+n/60]);
        h(n+42).LineWidth=10;
        c1=255;
        c2=0;
        c3=0;
        h(n+42).Color=[(c1-(n-1)*5.25)/255 c2/255 c3/255];
    end
    o=ones(1,length(K));
    p=ones(1,length(K));
    q=ones(1,length(K));
    o(1)=255;
    p(1)=255;
    q(1)=255;
    for n=2:length(K)
        o(n)=o(n-1)-7.419;
        p(n)=p(n-1)-7.419;
        q(n)=q(n-1)-4.6129;
    end
%     for n=1:length(K)        
%         %h(n+62)=plot([MIN_X*1.05+0.98*Kgradient/length(K)*(n-1) MIN_X*1.05+0.98*Kgradient/length(K)*n],[0 0]);
%         h(n+62)=plot([MIN_X*1.1+(Kgradient)*0.95-0.95*Kgradient/length(K)*(n-1) MIN_X*1.1+(Kgradient)*0.95-0.95*Kgradient/length(K)*n],[0 0]);
%         h(n+62).LineWidth=8;
%         h(n+62).Color=[o(n)/255 p(n)/255 q(n)/255];
%     end
% h=colorbar('southoutside')
% h.Limits = [Tmin Tmax];
frame = getframe(gcf);
hold off
close all
writeVideo(v,frame);
end
close(v)


end
