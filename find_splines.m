function [A B C D] = find_splines(X, Y, fpl, fpr)
    % compute coefficients in a(x-x0)^3+... form
    N = size(X,2) - 1;
    H = diff(X);
    alpha(1) = 3 * (Y(2) - Y(1)) / H(1) - 3 * fpl;
    alpha(N+1) = 3 * fpr - 3 * (Y(end) - Y(end-1)) / H(end);
    for i = 2:N
        alpha(i) = 3 * (Y(i+1) - Y(i)) / H(i) - 3 * (Y(i) - Y(i-1)) / H(i-1);
    end
    tridiag = zeros(N+1,N+1);
    diag_elts(N+1) = 2 * H(end);
    diag_elts(1) = 2 * H(1);
    diag_elts(2:N) = 2 * (H(1:N-1) + H(2:N));
    tridiag(1:N+2:(N+1)*(N+1)) = diag_elts; % diagonal elements
    tridiag(2:N+2:(N+1)*(N+1)) = H; % off-diagonal elements
    tridiag(N+2:N+2:(N+1)*(N+1)) = H;
    Bshift = (tridiag \ alpha')';
    for j = 1:N
        Cshift(j) = (Y(j+1) - Y(j)) / H(j) - H(j) * (Bshift(j+1) + 2 * Bshift(j)) / 3;
        Ashift(j) = (Bshift(j+1) - Bshift(j)) / (3 * H(j));
    end
    Dshift = Y(1:end-1);
    Bshift = Bshift(1:end-1); % trim off the last entry
    % Translate back to ax^3+bx^2+cx+d form
    A = Ashift;
    B = Bshift - 3 * Ashift .* X(1:end-1);
    C = 3 * Ashift .* X(1:end-1).^2 - 2 * Bshift .* X(1:end-1) + Cshift;
    D = Dshift - Ashift .* X(1:end-1).^3 + Bshift .* X(1:end-1).^2 - Cshift .* X(1:end-1);
end