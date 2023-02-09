function X = newt_vec(f, X0, TolFun, TolX, varargin)
    max_iter = 400; % maximum number of iterations
    status = 0; % 0 = stop within function tolerance, 1 = stop within step size tolerance, 2 = exceeded maximum number of iterations
    g = @(x) f(x, varargin{:});
    X = X0;
    fval = g(X);
    iter = 0;
    step = 0;
%     while norm(fval) > TolFun * (1 + norm(fval))
    while sum(fval.^2) > sqrt(TolFun)
%     while sum((g(X) - g(X + step')).^2) > TolFun * (1 + sum(fval.^2))
        iter = iter + 1;
        if iter > max_iter
            status = 2;
            break;
        end
%         step = g(X) * inv(jac(g,X))'; % g(X) is a row vector
        step = inv(jac(g,X)) * g(X); % g(X) is a column vector
%         if norm(step) < TolX * (1 + norm(X))
%             status = 1;
%             break;
%         end
        X = X - step'; % g(X) is a column vector but X is a row vector. why did I do this
        fval = g(X);
    end
    switch status
        case 0
            disp('Solver exited; function value within tolerance.');
%         case 1
%             disp('Solver exited; step size within tolerance.');
        case 2
            disp('Solver failed; exceeded maximum number of iterations.');
    end
end

function J = jac(f, X)
    N = size(X, 2);
    h = eps^(1/3);
    for i = 1:N
        % central differences
        J(:,i) = (f([X(1:i-1), X(i)+h, X(i+1:end)]) - f([X(1:i-1), X(i)-h, X(i+1:end)])) / (2*h);
    end
end