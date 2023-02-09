%% Compute stable states
% external force parameters
global Fext external_force_mode variance_factor
external_force_mode = "gaussian"; % one of "tip_only" or "gaussian"
Fext = 0.2;
variance_factor = 1/2;
initsize = 4;
endsize = 7;
results = cell(endsize-initsize+1, 1);

for logsize=initsize:endsize
    % set up the row in the results matrix for this iteration
    N = 2^logsize;
    % parameters and setup
    global Tol Rtol TolFun TolX Inc L P
    %dimensional parameters
    L = 1;%characteristic length
    P = 1;%characteristic pressure
    % Bulk Modulus
    k = 1;
    % Shear Modulus
    mu = 1;
    % Uniform material property
    K = k*ones(N,1);
    MU = mu*ones(N,1);
    %computational parameters
    Rtol = 1e-13;
    Tol =  1e-3;
    TolFun = 1e-16;
    TolX = 1e-16;
    Inc =  1;
    % Initialize cells
    R = 1;%the intrinsic radius, if a sphere is to be generated
    L0 = ones(N,1);
    R0 = ones(N,1);
    [ rv, adj ] = generate_sphere(R,N);
    %[ rv, adj ] = generate_arc(N);
    LUT = zeros(N,size(adj,1));
    ii = 1;
    for n = 1:N
        nv = find(adj(n,:));
        nv = nv(nv>n);
        for nn = nv
            LUT(ii,n) = 1; LUT(ii,nn) = -1;
            ii = ii + 1;
        end
    end
    for i=N:-1:1
       index = find(LUT(i,:));
       %R0(i) = min((rv(2,index(1))+rv(2,index(2)))/2,0.3*r0);
       L0(i) = sqrt((rv(1,index(1))-rv(1,index(2))).^2+(rv(2,index(1))-rv(2,index(2))).^2);
       R0(i) = (rv(2,index(1))+rv(2,index(2)))/2;
    end
    %get the bond to vertex map
    nbonds = sum(adj(:))/2;
    Tl = zeros(nbonds,1);% meridional tension
    Tr = zeros(nbonds,1);% circumferential tension
    ext_verts = find(sum(adj,1)==1);

    % find the stable configuration: first with just pressure
    [ Tl,Tr,rv,error ] = equilibrium(rv,LUT,L0,R0,K,MU,ext_verts,0);
    % then also with external force
    [ Tl,Tr,rv,error ] = equilibrium(rv,LUT,L0,R0,K,MU,ext_verts,1);
    
    % save the positions in the results matrix
    results(logsize-initsize+1) = {rv};
end

save(['erroranalysis_min_', num2str(initsize), '_max_', num2str(endsize), '_Fext_', num2str(Fext), '_var_', num2str(variance_factor), '.mat'], 'results');

%% Compute and plot positional error
results_trimmed{1} = results{1};
for stepsize=1:endsize-initsize
    results_trimmed{stepsize+1} = results{stepsize+1}(:,1:2^stepsize:end);
end

for logsize=initsize:endsize
    %pos_error{logsize-initsize+1} = results_trimmed{logsize-initsize+1}(:,:) - results_trimmed{1}(:,:);
    pos_error{logsize-initsize+1} = results_trimmed{logsize-initsize+1}(:,:) - results_trimmed{max(1,logsize-initsize)}(:,:);
    pos_error_norm{logsize-initsize+1} = vecnorm(pos_error{logsize-initsize+1}(:,:),2);
end

% plot(results_trimmed{1}(1,:), results_trimmed{1}(2,:), results_trimmed{2}(1,:), results_trimmed{2}(2,:), results_trimmed{3}(1,:), results_trimmed{3}(2,:), results_trimmed{4}(1,:), results_trimmed{4}(2,:)); pbaspect([1 1 1]);

for logsize = initsize+1:endsize
    overall_error_norm(1,logsize-initsize) = norm(pos_error_norm{logsize-initsize+1}, 1);
    overall_error_norm(2,logsize-initsize) = norm(pos_error_norm{logsize-initsize+1}, 2);
    overall_error_norm(3,logsize-initsize) = norm(pos_error_norm{logsize-initsize+1}, Inf);
end

plot([32, 64, 128], overall_error_norm(1,:), ...
    [32, 64, 128], overall_error_norm(2,:), ...
    [32, 64, 128], overall_error_norm(3,:)); ...
    set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');

saveas(gcf,['erroranalysis_position_min_', num2str(initsize), '_max_', num2str(endsize), '_Fext_', num2str(Fext), '_var_', num2str(variance_factor), '.png'])

%% Order of convergence for positional error
pos_order = log(norm(pos_error_norm{3}, Inf) / norm(pos_error_norm{4}, Inf)) / log(2);
fprintf("The computed order of convergence for positions is: %f.\n", pos_order);

%% Compute curvatures; compute and plot curvature error
for i=1:endsize-initsize+1
    [curvatures{i}(1,:) curvatures{i}(2,:)] = compute_curvatures(results_trimmed{i});
end

for logsize=initsize:endsize
    curv_ks_error{logsize-initsize+1} = curvatures{logsize-initsize+1}(1,:) - curvatures{max(1,logsize-initsize)}(1,:);
    curv_kphi_error{logsize-initsize+1} = curvatures{logsize-initsize+1}(2,:) - curvatures{max(1,logsize-initsize)}(2,:);
end

for logsize = initsize+1:endsize
    overall_curv_ks_error_norm(1,logsize-initsize) = norm(curv_ks_error{logsize-initsize+1}, 1);
    overall_curv_ks_error_norm(2,logsize-initsize) = norm(curv_ks_error{logsize-initsize+1}, 2);
    overall_curv_ks_error_norm(3,logsize-initsize) = norm(curv_ks_error{logsize-initsize+1}, Inf);
    overall_curv_kphi_error_norm(1,logsize-initsize) = norm(curv_kphi_error{logsize-initsize+1}, 1);
    overall_curv_kphi_error_norm(2,logsize-initsize) = norm(curv_kphi_error{logsize-initsize+1}, 2);
    overall_curv_kphi_error_norm(3,logsize-initsize) = norm(curv_kphi_error{logsize-initsize+1}, Inf);
end

% plot error for ks
plot([32, 64, 128], overall_curv_ks_error_norm(1,:), ...
    [32, 64, 128], overall_curv_ks_error_norm(2,:), ...
    [32, 64, 128], overall_curv_ks_error_norm(3,:)); ...
    set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
saveas(gcf,['erroranalysis_ks_min_', num2str(initsize), '_max_', num2str(endsize), '_Fext_', num2str(Fext), '_var_', num2str(variance_factor), '.png'])

% plot error for kphi
plot([32, 64, 128], overall_curv_kphi_error_norm(1,:), ...
    [32, 64, 128], overall_curv_kphi_error_norm(2,:), ...
    [32, 64, 128], overall_curv_kphi_error_norm(3,:)); ...
    set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
saveas(gcf,['erroranalysis_kphi_min_', num2str(initsize), '_max_', num2str(endsize), '_Fext_', num2str(Fext), '_var_', num2str(variance_factor), '.png'])

%% Order of convergence for curvatures
ks_order = log(norm(curv_ks_error{3}, Inf) / norm(curv_ks_error{4}, Inf)) / log(2);
fprintf("The computed order of convergence for ks is: %f.\n", ks_order);
kphi_order = log(norm(curv_kphi_error{3}, Inf) / norm(curv_kphi_error{4}, Inf)) / log(2);
fprintf("The computed order of convergence for kphi is: %f.\n", kphi_order);