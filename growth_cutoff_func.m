f0 = 2;
step = 0.0001;
invStep = int32(1/step);
X = 0:step:2;
F = zeros(1,2*invStep+1);
M = 1.487;
F = (1 - cutoff(X / M)) .* (cos(pi*X/M) + 1) / 2 + 1;
% for i = 0:step:2*invStep
%     F(i) = 1;
% end
Zprime = zeros(1,invStep+1);
Z = zeros(1,invStep+1);
Zprime(1) = 0;
Z(1) = 0;
for i = 1:1/step
    Zprime(i+1) = sqrt(F(i)^2 + 2 * (-1 + i * step) * F(i) * Fprime(i) + (-2 + i * step) * i * step * Fprime(i)^2);
    Z(i+1) = step * sum(Zprime(1:i));
%     Zprime(i+1) = sqrt(F(i)^2 - (Fprime(i) * sin(i * step) + F(i) * cos(i * step))^2);
%     Z(i+1) = step * Zprime(i) + Z(i);
end
plot(Z, F(1:invStep+1) .* sqrt(1 - ((0:step:1) - 1).^2) / sqrt(2));
% plot(Z / 1.182, F(1:invStep+1) .* sqrt(1 - ((0:step:1) - 1).^2) / 1.528);
hold on;
plot(X, sqrt(1 - (X - 1).^2));
daspect([1 1 1]);

function c = cutoff(x)
    if x <= 0
        c = 0;
    elseif x >= 1
        c = 1;
    else
        c = (1 + exp(1./x - 1./(1-x))).^-1;
    end
end