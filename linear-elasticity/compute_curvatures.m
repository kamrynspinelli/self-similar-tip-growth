function [ ks, kphi ] = compute_curvatures(rv)
tempVector = ones(2, 1);
ks = ones(1, size(rv, 2)-1);
kphi = ones(1, size(rv, 2)-1);
alpha = ones(1, size(rv, 2)-1);
dalpha = ones(1, size(rv, 2)-1);
ds = ones(1, size(rv, 2)-1);
rm = 0.5*(rv(2,1:end-1)+rv(2,2:end));
sigmass = ones(1, size(rv, 2)-1);
sigmaphiphi = ones(1, size(rv, 2)-1);
for i=1:(size(rv, 2)-1)
    tempVector(1, 1) = rv(1, i+1) - rv(1, i);
    tempVector(2, 1) = rv(2, i+1) - rv(2, i);
    ds(1, i) = sqrt(tempVector(1, 1)^2+tempVector(2, 1)^2);
    alpha(1, i) = atan(abs(tempVector(2, 1))/tempVector(1,1));
    kphi(1, i) = cos(alpha(1, i))/rm(i);
end
for i=1:(size(rv, 2)-2)
    dalpha(1, i) = alpha(1, i+1)-alpha(1, i);
end
dalpha(1, size(rv, 2)-1)=pi-2*alpha(1, size(rv, 2)-1);
for i=1:(size(rv, 2)-1)
    ks(1, i) = dalpha(1, i)/ds(1, i);
end
% for i=1:(size(rv, 2)-1)
%     sigmass(1, i) = 0.5/kphi(1, i);
%     sigmaphiphi(1, i) = (1-ks(1, i)*sigmass(1, i))/kphi(1, i);
% end