load('../cell-profiles/dumais-2004/9_2dF_cell1.mat'); % load the cell 1 data Dumais sent us
tipX = [678.2824  645.2651  612.3227  579.4825  545.8388];
tipY = [120.6311  124.3285  121.1008  116.6520  116.7097];
for i = size(Px,2):-5:size(Px,2)-20
    X = Px(1:max(find(Px(:,i))),i)'; % load the x and y coordinates
    Y = Py(1:max(find(Px(:,i))),i)'; % they trace the profile from bottom left, to tip, to top left
%     Xp = polyfit(floor(size(X,2))/6:5*floor(size(X,2))/6, X(floor(size(X,2))/6:5*floor(size(X,2))/6), 4);
%     Yp = polyfit(floor(size(X,2))/6:5*floor(size(X,2))/6, Y(floor(size(X,2))/6:5*floor(size(X,2))/6), 4);
%     Xpoly = @(t) Xp(1) * t.^4 + Xp(2) * t.^3 + Xp(3) * t.^2+ Xp(4) * t + Xp(5);
%     Ypoly = @(t) Yp(1) * t.^4 + Yp(2) * t.^3 + Yp(3) * t.^2+ Yp(4) * t + Yp(5);
%     hold on; 
%     plot(Xpoly(floor(size(X,2))/6:0.01:5*floor(size(X,2))/6), Ypoly(floor(size(X,2))/6:0.01:5*floor(size(X,2))/6));
%     plot(X, Y, 'x');
%     daspect([1 1 1]);
%     close all;
    
%     degree = 4
%     Xp = polyfit(Y(floor(size(X,2))/6:5*floor(size(X,2))/6), X(floor(size(X,2))/6:5*floor(size(X,2))/6), 4);
%     Xpoly = @(t) Xp(1) * t.^4 + Xp(2) * t.^3 + Xp(3) * t.^2+ Xp(4) * t + Xp(5);
%     degree = 6
%     Xp = polyfit(Y(floor(size(X,2))/6:5*floor(size(X,2))/6), X(floor(size(X,2))/6:5*floor(size(X,2))/6), 6);
%     Xpoly = @(t) Xp(1) * t.^6 + Xp(2) * t.^5 + Xp(3) * t.^4+ Xp(4) * t.^3 + Xp(5) * t.^2 + Xp(6) * t + Xp(7);
%     degree = 8
%     Xp = polyfit(Y(floor(size(X,2))/6:5*floor(size(X,2))/6), X(floor(size(X,2))/6:5*floor(size(X,2))/6), 8);
%     Xpoly = @(t) Xp(1) * t.^8 + Xp(2) * t.^7 + Xp(3) * t.^6 + Xp(4) * t.^5 + Xp(5) * t.^4 + Xp(6) * t.^3 + Xp(7) * t.^2 + Xp(8) * t + Xp(9);
%     degree = 10
    Xp = polyfit(Y(floor(size(X,2))/6:5*floor(size(X,2))/6), X(floor(size(X,2))/6:5*floor(size(X,2))/6), 10);
    Xpoly = @(t) Xp(1) * t.^10 + Xp(2) * t.^9 + Xp(3) * t.^8 + Xp(4) * t.^7 + Xp(5) * t.^6 + Xp(6) * t.^5 + Xp(7) * t.^4 + Xp(8) * t.^3 + Xp(9) * t.^2 + Xp(10) * t + Xp(11);
    hold on; 
    plot(Xpoly_all(Y(floor(size(X,2)/6)):0.01:Y(floor(5*size(X,2)/6))), Y(floor(size(X,2)/6)):0.01:Y(floor(5*size(X,2)/6)));
    plot(X, Y, 'x');
%     daspect([1 1 1]);
%     close all;
end