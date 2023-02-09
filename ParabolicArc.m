classdef ParabolicArc
    % ParabolicArc
    % Describes a parametric parabola fitting four points q1, p1, p2, q2, 
    % and its geometric properties.
    
    properties
        f % the parametric arc
        df % its velocity vector
        n % its outward normal vector
        tmin % f(tmin) = p1
        tmax % f(tmax) = p2
        tmid % f(tmid) = vertex of the parabola
        arclength % the arclength of f between tmin and tmax
        b % the midpoint of the segment p1p2
        L1 % the parabolic stretch factor fitting p1, p2, and q1
        L2 % the parabolic stretch factor fitting p1, p2, and q2
        L % the harmonic mean of L1 and L2
        ell % half the length of the segment p1p2
        beta % the signed angle from the basis vector z to the outward 
        %   normal of the segment p1p2
        alpha % the signed angle from the linear segment p1p2 to the basis 
        %   vector r
        % Tz % the z-coordinate of the vertex of the parabola
        % Tr % the r-coordinate of the vertex of the parabola
        K % bulk modulus
        mu % shear modulus
        vert % the vertex of the arc (midpoint of the patch by arclength for tip-adjacent patch)
        sigmaS % the meridional tension on the patch
        sigmaTheta % the circumfrential tension on the patch
    end
    
    methods
        function obj = ParabolicArc(q1, p1, p2, q2, K, mu, charS, charTheta, isTip, isDegenerate, tipPt, isUpper)
            % ParabolicArc: Construct an instance of this class
            %   Accepts 2-vectors p1, p2, q1, q2 describing the (z,r)
            %   coordinates of these four points, K and mu are the bulk and 
            %   shear modulus respectively, and charS and charTheta are the
            %   intrinsic arclength of the patch and charTheta is the 
            %   intrinsic r-coordinate of the patch midpoint. Returns a 
            %   ParabolicArc object which fits p1 and p2 and compromises 
            %   between fitting q1 and q2.
            if nargin == 0 % if no arguments were supplied,
                obj.f = @(t) 0;
                obj.df = @(t) 0;
                obj.n = @(t) 0;
                obj.tmin = 0;
                obj.tmax = 0;
                obj.tmid = 0;
                obj.arclength = 0;
                obj.b = [0; 0];
                obj.L1 = 0;
                obj.L2 = 0;
                obj.L = 0;
                obj.ell = 0;
                obj.beta = 0;
                obj.alpha = 0;
                obj.K = 0;
                obj.mu = 0;
                obj.vert = 0;
                obj.sigmaS = 0;
                obj.sigmaTheta = 0;
                return; % return a dummy ParabolicArc object
            end         
            obj.b = (p1 + p2) / 2; 
            obj.beta = atan((p2(2) - p1(2)) / (p2(1) - p1(1))) + pi/2;
            obj.alpha = pi - obj.beta;
            obj.ell = norm(p1 - obj.b, 2);
            if isTip % if this is a tip patch
                if isUpper % if it's above the z-axis
                    obj.tmax = 0;
                    obj.tmin = -obj.ell;
                    % obj.tmid to be assigned later
                else
                    obj.tmax = obj.ell;
                    obj.tmin = 0;
                    % obj.tmid to be assigned later
                end
            else
                obj.tmax = obj.ell;
                obj.tmin = -obj.ell;
                obj.tmid = 0;
            end
            % p1_rot = [-obj.ell, 0];
            % p2_rot = [obj.ell, 0];
            q1_rot = [sin(obj.beta) * (q1(1) - obj.b(1)) - cos(obj.beta) * (q1(2) - obj.b(2)), ...
                cos(obj.beta) * (q1(1) - obj.b(1)) + sin(obj.beta) * (q1(2) - obj.b(2))];
            q2_rot = [sin(obj.beta) * (q2(1) - obj.b(1)) - cos(obj.beta) * (q2(2) - obj.b(2)), ...
                cos(obj.beta) * (q2(1) - obj.b(1)) + sin(obj.beta) * (q2(2) - obj.b(2))];
            if isDegenerate
                obj.L = Inf;
            elseif isTip
                obj.L = p1(2) / ParabolicArc.signed_sqrt(tipPt(1) - p1(1));
            else
                obj.L1 = ParabolicArc.signed_sqrt((obj.ell^2 - q1_rot(1)^2) / q1_rot(2));
                obj.L2 = ParabolicArc.signed_sqrt((obj.ell^2 - q2_rot(1)^2) / q2_rot(2));
                obj.L = (2 * obj.L1 * obj.L2) / (obj.L1 + obj.L2);
            end
            if p2(1) < p1(1) % if p2 is to the left of p1
                obj.f = @(t) [sin(obj.beta), cos(obj.beta); -cos(obj.beta), sin(obj.beta)] ...
                    * [-t ; (-t.^2 + obj.ell^2) / (sign(obj.L) * obj.L^2)] ...
                    + obj.b; % use -t instead of t to parametrize in the right direction
                obj.df = @(t) [sin(obj.beta), cos(obj.beta); -cos(obj.beta), sin(obj.beta)] ...
                    * [-1 ; (-2 * t) / (sign(obj.L) * obj.L^2)];
%                 obj.n = @(t) [0, -1; 1, 0] ... % rotate by pi/2 counterclockwise
%                     * [sin(obj.beta), cos(obj.beta); -cos(obj.beta), sin(obj.beta)] ...
%                     * [-1 ; (-2 * t) / (sign(obj.L) * obj.L^2)];
                obj.n = @(t) [-cos(obj.beta) + sin(obj.beta) * 2 * t / (sign(obj.L) * obj.L^2); ...
                    -sin(obj.beta) - cos(obj.beta) * 2 * t / (sign(obj.L) * obj.L^2)];
            else
                obj.f = @(t) [sin(obj.beta), cos(obj.beta); -cos(obj.beta), sin(obj.beta)] ...
                    * [t ; (-t.^2 + obj.ell^2) / (sign(obj.L) * obj.L^2)] ...
                    + obj.b;
                obj.df = @(t) [sin(obj.beta), cos(obj.beta); -cos(obj.beta), sin(obj.beta)] ...
                    * [1 ; (-2 * t) / (sign(obj.L) * obj.L^2)];
%                 obj.n = @(t) [0, -1; 1, 0] ... % rotate by pi/2 counterclockwise
%                     * [sin(obj.beta), cos(obj.beta); -cos(obj.beta), sin(obj.beta)] ...
%                     * [1 ; (-2 * t) / (sign(obj.L) * obj.L^2)];
                obj.n = @(t) [cos(obj.beta) + sin(obj.beta) * 2 * t / (sign(obj.L) * obj.L^2); ...
                    sin(obj.beta) - cos(obj.beta) * 2 * t / (sign(obj.L) * obj.L^2)];
            end
            if abs(obj.L) == Inf % covers degenerate and flat cases
                obj.arclength = 2 * obj.ell;
            elseif isTip % if this is a tip patch, take the arclength from one side of the vertex
                obj.arclength = 1/4 * (2 * obj.ell * sqrt(4 * obj.ell^2 / obj.L^4 + 1) ...
                    + obj.L^2 * asinh(2 * obj.ell / obj.L^2));
            else
                obj.arclength = 1/2 * (2 * obj.tmax * sqrt(4 * obj.tmax^2 / obj.L^4 + 1) ...
                    + obj.L^2 * asinh(2 * obj.tmax / obj.L^2));
            end
            if isTip
                absMidpt = midpoint_by_arclength(@(t) 2*t / obj.L^2, 0, obj.ell); % the absolute value of t giving a midpoint of the patch by arclength
                if isUpper
                    obj.tmid = -absMidpt;
                else
                    obj.tmid = absMidpt;
                end
            end
            obj.K = K;
            obj.mu = mu;
            obj.vert = obj.f(obj.tmid); % the vertex of the parabola
            Ls = obj.arclength / charS; % the arclength ratio ds/ds^0
            LsInv = charS / obj.arclength; % the arclength ratio ds^0/ds
            Ltheta = obj.vert(2) / charTheta; % the circumferential ratio r/r_0
            LthetaInv = charTheta / obj.vert(2); % the circumferential ratio r_0/r
            obj.sigmaS = 1/2 * mu * (LthetaInv^2 - LsInv^2) ...
                + K * (Ls * Ltheta - 1);
            obj.sigmaTheta = 1/2 * mu * (LsInv^2 - LthetaInv^2) ...
                + K * (Ls * Ltheta - 1);
%             obj.sigmaS = 1/2 * mu * (Ls^2 - Ltheta^2) / (Ltheta^2 * Ls^2) + K * (Ls * Ltheta - 1);
%             obj.sigmaTheta = 1/2 * mu * (Ltheta^2 - Ls^2) / (Ltheta^2 * Ls^2) + K * (Ls * Ltheta - 1);
        end
    end
    
    methods (Static = true)
        function s = signed_sqrt(a)
            s = sign(a) * sqrt(abs(a));
        end
        
        function arcs = all_arcs(points, K, mu, charS, charTheta)
           % all_arcs
           %    Accepts a 2xN matrix points whose colums are marker points, 
           %    in order; N-1-vectors K and mu of bulk and shear moduli for 
           %    each patch, and N-1-vectors charS and charTheta containing 
           %    the intrinsic arclengths and circumferential radii of each
           %    patch. Returns an array of the N+1 ParabolicArc objects
           %    describing the patches between the points and the two 
           %    exterior patches adjacent to the start and end points, 
           %    assuming that the points all lie in the first quadrant with 
           %    the first point on the r-axis, and that the membrane is 
           %    symmetric with respect to the z- and r-axes.
           N = size(points, 2); % the number of points
           ext_points = zeros(2, N+4); % preallocate extended points matrix
           % fill the first two spots with the first two points, reflected
           ext_points(1, 1:2) = -fliplr(points(1, 2:3));
           ext_points(2, 1:2) = fliplr(points(2, 2:3));
           ext_points(:, 3:N+2) = points;
           ext_points(1, N+3:N+4) = fliplr(points(1, N-2:N-1));
           ext_points(2, N+3:N+4) = -fliplr(points(2, N-2:N-1));
           ext_K = zeros(1, N+1); % preallocate extended patch parameter array
           ext_K(1) = K(1); % fill in
           ext_K(2:N) = K(:);
           ext_K(N+1) = K(N-1);
           ext_mu = zeros(1, N+1);
           ext_mu(1) = mu(1);
           ext_mu(2:N) = mu(:);
           ext_mu(N+1) = mu(N-1);
           ext_charS = zeros(1, N+1);
           ext_charS(1) = charS(1);
           ext_charS(2:N) = charS(:);
           ext_charS(N+1) = charS(N-1);
           ext_charTheta = zeros(1, N+1);
           ext_charTheta(1) = charTheta(1);
           ext_charTheta(2:N) = charTheta(:);
           ext_charTheta(N+1) = -charTheta(N-1); % need to negate this one
           % preallocate object array by assigning last element
           arcs(N+1) = ParabolicArc(ext_points(:,N), ext_points(:,N+1), ...
               ext_points(:,N+3), ext_points(:,N+4), ext_K(N+1), ext_mu(N+1), ...
               ext_charS(N+1), ext_charTheta(N+1), true, false, ext_points(:,N+2), false);
           arcs(N) = ParabolicArc(ext_points(:,N), ext_points(:,N+1), ...
               ext_points(:,N+3), ext_points(:,N+4), ext_K(N), ext_mu(N), ...
               ext_charS(N), ext_charTheta(N), true, false, ext_points(:,N+2), true);
           % fill in the rest
           for i = 1:N-1
               arcs(i) = ParabolicArc(ext_points(:,i), ext_points(:,i+1), ...
                   ext_points(:,i+2), ext_points(:,i+3), ext_K(i), ext_mu(i), ...
                   ext_charS(i), ext_charTheta(i), false, false);
           end
        end
        
        function arcs = all_degenerate_arcs(points, K, mu, charS, charTheta)
           % all_arcs
           %    Accepts a 2xN matrix points whose colums are marker points, 
           %    in order; N-1-vectors K and mu of bulk and shear moduli for 
           %    each patch, and N-1-vectors charS and charTheta containing 
           %    the intrinsic arclengths and circumferential radii of each
           %    patch. Returns an array of the N+1 degenerate 
           %    (\lambda = \infty)ParabolicArc objects describing the 
           %    patches between the points and the two exterior patches 
           %    adjacent to the start and end points, assuming that the 
           %    points all lie in the first quadrant with the first point 
           %    on the r-axis, and that the membrane is symmetric with 
           %    respect to the z- and r-axes.
           N = size(points, 2); % the number of points
           ext_points = zeros(2, N+4); % preallocate extended points matrix
           % fill the first two spots with the first two points, reflected
           ext_points(1, 1:2) = -fliplr(points(1, 2:3));
           ext_points(2, 1:2) = fliplr(points(2, 2:3));
           ext_points(:, 3:N+2) = points;
           ext_points(1, N+3:N+4) = fliplr(points(1, N-2:N-1));
           ext_points(2, N+3:N+4) = -fliplr(points(2, N-2:N-1));
           ext_K = zeros(1, N+1); % preallocate extended patch parameter array
           ext_K(1) = K(1); % fill in
           ext_K(2:N) = K(:);
           ext_K(N+1) = K(N-1);
           ext_mu = zeros(1, N+1);
           ext_mu(1) = mu(1);
           ext_mu(2:N) = mu(:);
           ext_mu(N+1) = mu(N-1);
           ext_charS = zeros(1, N+1);
           ext_charS(1) = charS(1);
           ext_charS(2:N) = charS(:);
           ext_charS(N+1) = charS(N-1);
           ext_charTheta = zeros(1, N+1);
           ext_charTheta(1) = charTheta(1);
           ext_charTheta(2:N) = charTheta(:);
           ext_charTheta(N+1) = -charTheta(N-1); % need to negate this one
           % preallocate object array by assigning last element
           arcs(N+1) = ParabolicArc(ext_points(:,N+1), ext_points(:,N+2), ...
               ext_points(:,N+3), ext_points(:,N+4), ext_K(N+1), ext_mu(N+1), ...
               ext_charS(N+1), ext_charTheta(N+1), false, true);
           % fill in the rest
           for i = 1:N
               arcs(i) = ParabolicArc(ext_points(:,i), ext_points(:,i+1), ...
                   ext_points(:,i+2), ext_points(:,i+3), ext_K(i), ext_mu(i), ...
                   ext_charS(i), ext_charTheta(i), false, true);
           end
           for i = 1:N+1
               arcs(i).L = Inf;
           end
        end
        
        function p = plot_arcs(points, heatVals, heatMin, heatMax)
            % plot_arcs
            %   Accepts a 2xN matrix points whose columns are marker
            %   points. Plots all the parabolic arcs for the patches
            %   between these marker points
            if nargin >= 2
                cmap = jet;
            end
            N = size(points, 2); % the number of points
            arcs = ParabolicArc.all_arcs(points, ones(1, N-1), ones(1, N-1), ...
                ones(1, N-1), ones(1, N-1)); % compute all the arcs
            arcs = arcs(2:end-1); % remove the arcs that live outside the first quadrant
            hold on;
            for i = 1:N-1
                thisArc = arcs(i);
                plotData = thisArc.f(linspace(thisArc.tmin, thisArc.tmax, 20));
                p = plot(plotData(1,:), plotData(2,:), 'LineWidth', 1.0);
                if nargin >= 2
                    if heatMin == heatMax
                        colorVal = 1;
                    else
                        colorVal = min(max(round((heatVals(i) - heatMin) / (heatMax - heatMin) * 256), 1), 256);
                    end
                    p.Color = cmap(colorVal, :);
                end
            end
            daspect([1 1 1]);
            if nargin >= 2
                colormap jet; colorbar; caxis([heatMin heatMax]);
            end
            p = gcf;
        end
        
        function p = plot_degenerate_arcs(points)
            % plot_degenerate_arcs
            %   Accepts a 2xN matrix points whose columns are marker
            %   points. Plots all the degenerate arcs for the patches
            %   between these marker points
            N = size(points, 2); % the number of points
            arcs = ParabolicArc.all_degenerate_arcs(points, ones(1, N-1), ones(1, N-1), ...
                ones(1, N-1), ones(1, N-1)); % compute all the arcs
            arcs = arcs(2:end-1); % remove the arcs that live outside the first quadrant
            hold on;
            for i = 1:N-1
                thisArc = arcs(i);
                plotData = thisArc.f(linspace(thisArc.tmin, thisArc.tmax, 20));
                plot(plotData(1,:), plotData(2,:), 'LineWidth', 1.0);
            end
            daspect([1 1 1]);
            p = gcf;
        end
        
        function p = plot_arcs_with_displacement(approx, exact)
            % plot_arcs
            %   Accepts a 2xN matrix approx whose columns are marker
            %   points of an approximate profile and a 2xN matrix exact 
            %   whose columns are marker points of an exact profile. Plots 
            %   all the parabolic arcs for the exact and approximate
            %   profiles, and the displacement vectors from the exact
            %   marker points to the approximate ones.
            N = size(approx, 2);
            displacements = approx - exact;
            p = ParabolicArc.plot_arcs(approx);
            hold on;
            exactArcs = ParabolicArc.all_arcs(exact, ones(1, N-1), ones(1, N-1), ...
                ones(1, N-1), ones(1, N-1)); % compute all the arcs for the exact profile
            exactArcs = exactArcs(2:end-1); % remove the arcs that live outside the first quadrant
            hold on;
            for i = 1:N-1
                thisArc = exactArcs(i);
                plotData = thisArc.f(linspace(thisArc.tmin, thisArc.tmax, 20));
                plot(plotData(1,:), plotData(2,:), 'LineWidth', 1.0);
            end
            quiver(exact(1,:), exact(2,:), displacements(1,:), displacements(2,:), 0);% plot the displacements
            daspect([1 1 1]);
            p = gcf;
        end
        
        function dot = Tvec2D(point, lPatch, rPatch)
            % Tvec2D
            %   Accepts a column 2-vector point representing the 
            %   coordinates of a point and two ParabolicArc objects lPatch 
            %   and rPatch which are the patches preceding and following 
            %   the marker point. Returns the component of force on the
            %   marker point corresponding to S_1 in the writeup.
            lTvec = lPatch.df(lPatch.tmid); % a tangent vector to the left-side patch at its midpoint
            lTvec = lTvec / norm(lTvec); % make it a unit vector
            rTvec = rPatch.df(rPatch.tmid); % a tangent vector to the left-side patch at its midpoint
            rTvec = rTvec / norm(rTvec); % make it a unit vector
            lSigmaS = lPatch.sigmaS;
            rSigmaS = rPatch.sigmaS;
            dot = rSigmaS * rTvec - lSigmaS * lTvec;
        end
        
        function dot = Tvec3DrR(point, lPatch, rPatch)
            % Tvec3DrR
            %   Accepts a column 2-vector point representing the 
            %   coordinates of a point and two ParabolicArc objects lPatch 
            %   and rPatch which are the patches preceding and following 
            %   the marker point. Returns the component of force on the 
            %   marker point corresponding to S_2 in the writeup.
            lSigmaS = lPatch.sigmaS; % needed patch parameters
            rSigmaS = rPatch.sigmaS;
            lSigmaTheta = lPatch.sigmaTheta;
            rSigmaTheta = rPatch.sigmaTheta;
            lAlpha = lPatch.alpha;
            rAlpha = rPatch.alpha;
            lEll = lPatch.ell;
            rEll = rPatch.ell;
            lLambda = lPatch.L;
            rLambda = rPatch.L;
            lMdpt = lPatch.b;
            rMdpt = rPatch.b;
            lTmax = lPatch.tmax;
            lTmid = lPatch.tmid;
            rTmin = rPatch.tmin;
            rTmid = rPatch.tmid;
            % integrands for the left- and right-side integrals
            lIntegrand = @(t) (lSigmaS - lSigmaTheta) ...
                / (t * cos(lAlpha) + (-t^2 + lEll^2) / (sign(lLambda) * lLambda^2) * sin(lAlpha) + lMdpt(2)) ...
                * sqrt(1 + 4*t^2 / lLambda^4);
            rIntegrand = @(t) (rSigmaS - rSigmaTheta) ...
                / (t * cos(rAlpha) + (-t^2 + rEll^2) / (sign(rLambda) * rLambda^2) * sin(rAlpha) + rMdpt(2)) ...
                * sqrt(1 + 4*t^2 / rLambda^4);
            % compute the force
            dot = [0; midpt(lIntegrand, lTmid, lTmax, 16) ...
                + midpt(rIntegrand, rTmin, rTmid, 16)];
        end
        
        function dot = Tvec3DnR(point, lPatch, rPatch)
            % Tvec3DrR
            %   Accepts a column 2-vector point representing the 
            %   coordinates of a point and two ParabolicArc objects lPatch 
            %   and rPatch which are the patches preceding and following 
            %   the marker point. Returns the component of force on the 
            %   marker point corresponding to S_3 in the writeup.
            lSigmaS = lPatch.sigmaS; % needed patch parameters
            rSigmaS = rPatch.sigmaS;
            lAlpha = lPatch.alpha;
            rAlpha = rPatch.alpha;
            lEll = lPatch.ell;
            rEll = rPatch.ell;
            lLambda = lPatch.L;
            rLambda = rPatch.L;
            lMdpt = lPatch.b;
            rMdpt = rPatch.b;
            lNvec = lPatch.n;
            rNvec = rPatch.n;
            lTmax = lPatch.tmax;
            lTmid = lPatch.tmid;
            rTmin = rPatch.tmin;
            rTmid = rPatch.tmid;
            lMidTvec = lPatch.df(0);
            rMidTvec = rPatch.df(0);
%             if abs(point(2)) < 0.00001 % if the marker point is the tip,
%                 dot = [-lTmax; 0]; % assume that r = 2 \sigma_s \sin\alpha
%                 return; % TODO: double-check that this is reasonable
%             end
            % integrands for the left- and right-side integrals
            if lMidTvec(1) < 0 % if the left patch goes right-to-left
                lIntegrand = @(t) -lSigmaS * (-sin(lAlpha) + 2*t * cos(lAlpha) / (sign(lLambda) * lLambda^2)) ...
                    / (-t * cos(lAlpha) + (-t^2 + lEll^2) / (sign(lLambda) * lLambda^2) * sin(lAlpha) + lMdpt(2)) ...
                    * lNvec(t) / norm(lNvec(t)); % need to negate the sin(alpha) computation
            else
                lIntegrand = @(t) -lSigmaS * (sin(lAlpha) + 2*t * cos(lAlpha) / (sign(lLambda) * lLambda^2)) ...
                    / (t * cos(lAlpha) + (-t^2 + lEll^2) / (sign(lLambda) * lLambda^2) * sin(lAlpha) + lMdpt(2)) ...
                    * lNvec(t) / norm(lNvec(t));
            end
            if rMidTvec(1) < 0  % if the right patch goes right-to-left
                rIntegrand = @(t) -rSigmaS * (-sin(rAlpha) + 2 * t * cos(rAlpha) / (sign(rLambda) * rLambda^2)) ...
                    / (-t * cos(rAlpha) + (-t^2 + rEll^2) / (sign(rLambda) * rLambda^2) * sin(rAlpha) + rMdpt(2)) ...
                    * rNvec(t) / norm(rNvec(t)); % need to negate the sin(alpha) computation
            else
                rIntegrand = @(t) -rSigmaS * (sin(rAlpha) + 2*t * cos(rAlpha) / (sign(rLambda) * rLambda^2)) ...
                    / (t * cos(rAlpha) + (-t^2 + rEll^2) / (sign(rLambda) * rLambda^2) * sin(rAlpha) + rMdpt(2)) ...
                    * rNvec(t) / norm(rNvec(t));
            end
            % compute the force
            dot = midpt(lIntegrand, lTmid, lTmax, 16) ...
                + midpt(rIntegrand, rTmin, rTmid, 16);
        end
        
        function dot = Pvec(point, lPatch, rPatch, P)
            % Pvec
            %   Accepts a column 2-vector point representing the 
            %   coordinates of a point, two ParabolicArc objects lPatch 
            %   and rPatch which are the patches preceding and following 
            %   the marker point, and a scalar pressure P. Returns the
            %   component of force on the marker point corresponding to
            %   S_4 in the writeup.
            lNvec = lPatch.n; % the function handle for the normal to the left-side patch
            rNvec = rPatch.n; % the function handle for the normal to the left-side patch
            lLambda = lPatch.L; % the parabolic stretch factor on the left-side patch
            rLambda = rPatch.L; % the parabolic stretch factor on the left-side patch
            lTmax = lPatch.tmax;
            lTmid = lPatch.tmid;
            rTmin = rPatch.tmin;
            rTmid = rPatch.tmid;
            lIntgrnd = @(t) lNvec(t) / norm(lNvec(t)) ...
                * sqrt(1 + 4*t^2 / (lLambda^4)); % the integrand for the integral over the left half-patch
            rIntgrnd = @(t) rNvec(t) / norm(rNvec(t)) ...
                * sqrt(1 + 4*t^2 / (rLambda^4)); % the integrand for the integral over the right half-patch
            lInt = midpt(lIntgrnd, lTmid, lTmax, 16); % the numerical integration over the left half-patch
            rInt = midpt(rIntgrnd, rTmin, rTmid, 16); % the numerical integration over the left half-patch
            dot = P * (lInt + rInt);
        end
    end
end

