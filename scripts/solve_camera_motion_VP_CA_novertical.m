%% This code solves the minimization problem 

function [xopt, vpt, vp_asstn] = solve_camera_motion_VP_CA_novertical(linesRS, parm)

    poly_N = parm.poly_N;
     poly_funcX = @(coeff,t) coeff(1)*t.^2+coeff(2)*t; 
    poly_funcY = @(coeff,t) coeff(1)*t.^2+coeff(2)*t; 
    poly_funcZ = @(coeff,t) coeff(1)*t.^2+coeff(2)*t;
    h = parm.h; 
    invK = inv(parm.K);
    wt = compute_weights(linesRS);
    no_lines = size(linesRS, 1)/2; 

    tx = @(t) [t.^2, t]; 
    ty = @(t) [t.^2, t]; 
    tz = @(t) [t.^2, t]; 

    optmincon = optimoptions(@fmincon, 'MaxIter', 500, 'TolFun', 1.0000e-5); 
    optmincon = optimoptions(optmincon,'GradObj','off','GradConstr','on', 'Algorithm','interior-point');
    optmincon = optimoptions(optmincon, 'TolX', 1.0000e-5, 'Display',  'iter', 'Diagnostics', 'on'); 
    optmincon = optimoptions(optmincon, 'MaxFunEvals', 10000);%, 'Hessian', 'user-supplied', 'HessFcn', @hessinterior);%, 'TolX', 0, 'TolFun', 0);

    p = 3*(poly_N-1); 
    pr = [0.1*ones(1, floor((p-1)/(poly_N-1))); 0.05*ones(1, floor((p-1)/(poly_N-1)))];
    pr = 0*[0.1, 0.05, 0.1, 0.05, 0.1, 0.05];
    pr = pr(:); 
    
%     x0 = 0.02*randn(1, p); 
    x0 = 0.02*zeros(1, p); 
    x0 = [x0, 0, 0, 0]; 
    options.ub = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 0.25]';  % Lower bound on the variables.
    options.lb = - options.ub;  % Upper bound on the variables.

    hPixel = invK*[linesRS, ones(size(linesRS, 1), 1)]'; 
    
%     hPixelO = invK*[linesPts, ones(size(linesPts, 1), 1)]'; 
%     [xopt] = fmincon(@objective_fmincon, x0, A, Ab, [], [], options.lb, options.ub, [], optmincon);%@constraints_fmincon, optmincon);
    xoptm = fmincon(@objective_fmincon, x0, [], [], [], [], options.lb, options.ub, [], optmincon);
 
%     options = optimset('Display','iter', 'MaxFunEvals', 10000, 'MaxIter', 1500, 'TolFun',1e-12, 'TolX', 1.0000e-15);  
%     xoptm = fminsearch(@objective_fmincon,x0, options); 
    xopt.rx = xoptm(1:poly_N-1)'; 
    xopt.ry = xoptm(poly_N:2*(poly_N-1))'; 
    xopt.rz = xoptm(2*poly_N-1:3*(poly_N-1))'; 
     
    [vpt1, vpt2, vpt3] = cayley_rotation(xoptm(end-2:end));  
    vpt = [vpt1; vpt2; vpt3]';
    vp_asstn = vp_association(xoptm); 
    
% ----------------------------------------------------------------------
    function [f, g] = objective_fmincon (x)

        [vp1, vp2, vp3] = cayley_rotation(x(end-2:end));
        vp = [vp1; vp2; vp3]';
        
        i = 1:no_lines; 
        t = linesRS(2*i-1, 2)/h; 
        rxyz = [poly_funcX(x(1:poly_N-1), t), poly_funcY(x(poly_N:2*(poly_N-1)), t), poly_funcZ(x(2*poly_N-1:3*(poly_N-1)), t);];% poly_func(x(2*poly_N+1:3*poly_N), t)]; 
        [R1, R2, R3] = cayley_rotation(rxyz); 
        a1  = sum(R1'.*hPixel(:, 2*i-1));
        a2  = sum(R2'.*hPixel(:, 2*i-1));
        a3  = sum(R3'.*hPixel(:, 2*i-1));

        t = linesRS(2*i, 2)/h; 
        rxyz = [poly_funcX(x(1:poly_N-1), t), poly_funcY(x(poly_N:2*(poly_N-1)), t), poly_funcZ(x(2*poly_N-1:3*(poly_N-1)), t)];% poly_func(x(2*poly_N+1:3*poly_N), t)]; 
        [R1, R2, R3] = cayley_rotation(rxyz); 
        b1  = sum(R1'.*hPixel(:, 2*i));
        b2  = sum(R2'.*hPixel(:, 2*i));
        b3  = sum(R3'.*hPixel(:, 2*i));

        mHat = [a2.*b3-a3.*b2; a3.*b1-a1.*b3; a1.*b2-a2.*b1]; 
        mHat = mHat./repmat(sqrt(sum(mHat.^2)), [3, 1]); 

        vl = abs(vp*mHat); 
        fx = vl(2, :);
        id = vl(2, :) > parm.th;  
        [fx(id), ivd] = min(vl(:, id)); 
        f = huber(fx, parm.th)*wt; 

        ll = [((x(1:(poly_N-1):end-3))), ((x(2:(poly_N-1):end-3)))]; 
        f = sum(f)+abs(ll)*pr; % dont need any prior on x3
        g = []; 
    end

   
    function hx = huber(xx, th)
         ax = abs(xx);
         hx = 0.5*xx.^2;
         id = ax > th;  
         hx(id) = th*(ax(id) - 0.5*th);
    end

    function hx = tukey(xx, th)
         ax = abs(xx);
         hx = (th^2/6)*(1-(1-(xx./th).^2).^3);
         id = ax > th;  
         hx(id) = th^2/6;
    end

    function [R1, R2, R3] = cayley_rotation(rxyz)
            rx=rxyz(:, 1);
            ry=rxyz(:, 2); 
            rz=rxyz(:, 3);
            
            rxx=rx.*rx;
            ryy=ry.*ry;
            rzz=rz.*rz;
            
            rxy = rx.*ry; 
            rxz = rx.*rz; 
            ryz = ry.*rz; 
            
            Z = repmat((1+rxx+ryy+rzz).^(-1), [1, 3]); 
            R1 = Z.*[1+rxx-ryy-rzz, 2*(rxy+rz), 2*(rxz-ry)];
            R2 = Z.*[2*(rxy-rz), 1-rxx+ryy-rzz, 2*(ryz+rx)];
            R3 = Z.*[2*(rxz+ry), 2*(ryz-rx), 1-rxx-ryy+rzz;]; 
                
    end

    function vp_asstn = vp_association(x)
        [vp1, vp2, vp3] = cayley_rotation(x(end-2:end));
        vp = [vp1; vp2; vp3];
        vp_asstn = zeros(1, no_lines); 
        for i = 1:no_lines
            t = linesRS(2*i-1, 2)/h; 
            rxyz = [poly_funcX(x(1:poly_N-1), t); poly_funcY(x(poly_N:2*(poly_N-1)), t);poly_funcZ(x(2*poly_N-1:3*(poly_N-1)), t);];% poly_func(x(2*poly_N+1:3*poly_N), t)]; 
            [r1, r2, r3] = cayley_rotation(rxyz'); 
            r = [r1; r2; r3]; 
            a  = r*hPixel(:, 2*i-1);
            
            t = linesRS(2*i, 2)/h; 
            rxyz = [poly_funcX(x(1:poly_N-1), t); poly_funcY(x(poly_N:2*(poly_N-1)), t);poly_funcZ(x(2*poly_N-1:3*(poly_N-1)), t);];% poly_func(x(2*poly_N+1:3*poly_N), t)]; 
            [r1, r2, r3] = cayley_rotation(rxyz'); 
            r = [r1; r2; r3]; 
            b  = r*hPixel(:, 2*i);
             
            mHat = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]'; 
            mHat = mHat./norm(mHat); 
            
            vl = abs(vp'*mHat); 
            fx = vl(2, :);
            [mn, id] =  min(vl); 
            if fx < parm.th;
                vp_asstn(i) = 2;
            elseif mn < parm.th
                vp_asstn(i) = id; 
            end
        end
    end
end
  

function wt = compute_weights(linesRS)

    wt = zeros(size(linesRS, 1)/2, 1); 
    for i = 1:numel(wt)
        wt(i) = norm(linesRS(2*i-1, :) - linesRS(2*i, :)); 
    end
    wt = wt/max(wt); 
end

