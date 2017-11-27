function initial_shape = AFTER_GRAVITY(curvature,...
                                          torsion,...
                                          L,...
                                          gravity_on,...
                                          n_t,...
                                          n_f,...
                                          r_cyl)

    
    % Make structure to pass around parameters more easily
    g   = 9.80665;
    Eb  = 7e5;
    WpL = 0.0441;    % for 60, 88 actuator

    dib = 0.0095; % inner dia
    dob = 0.013;  % outer dia

    Ib   = pi/4*((dob/2)^4-(dib/2)^4);  % second moment wrt central axis
    W_b  = WpL*L;                       % weight of the actuator
    eg   = [0;0; gravity_on];           % global direction of gravity
    EC_F = 1e-3*g*eg;                   % End cap force (in global direction)
    rhoA = W_b/L;                       % mass per unit length
    C    = [1; 1; 2/3] * Eb*Ib;
    u_p  = [curvature                   % precurved tube u's
            0
            torsion];

    params = struct('u_p'          , u_p,...
                    'g'            , g,...
                    'C'            , C,...
                    'eg'           , eg,...
                    'rhoA'         , rhoA,...
                    'EC_F'         , EC_F,...
                    'R'            , eye(3),...
                    'gravity_flag' , gravity_on,...
                    'Eb'           , Eb,...
                    'L'            , L,...
                    'WpL'          , WpL,...
                    'n_t'          , n_t,...
                    'cyl_center'   , [0;0;0],... % to be optimized
                    'initial_shape', [],...
                    'r_cyl'        , r_cyl);

    % options for fmincon
%     options = optimoptions('fmincon'                ,...
%                            'Display'                ,'iter',...
%                            'Algorithm'              ,'sqp',... % use active-set or sqp only, interiror point doesn't satisfy the constraints
%                            'MaxFunctionEvaluations' ,10e3,...
%                            'MaxIter'                ,180,...
%                            'InitBarrierParam'       , 1e10,...
%                            'TolFun'                 , 1e-6,...
%                            'TolCon'                 , 1e-6,...
%                            'StepTolerance'          , 1e-10,...
%                            'UseParallel'            ,true);%,...
%     
    % Compute initial shape
    x0 = 0*ones(n_f,1);
    start_guess = [zeros(n_t,3)   ones(n_t,1)           zeros(n_t,3) ...
                   ones(n_t,1)    zeros(n_t,3)           ones(n_t,1) ...
                   zeros(n_t,2)   torsion*ones(n_t,1)   zeros(n_t,3)];

    params.initial_shape = start_guess;
    initial_shape = final_shape(x0, params);
    params.initial_shape = initial_shape;
    
    
    % Make initial plot
    clf, hold on
    plot3(params.initial_shape(:,1) ,...
          params.initial_shape(:,2),...
          params.initial_shape(:,3), ...
          'b');
    
   

    axis equal
    grid on
% clear all
    
end





function [sxint,...
          dist,...
          ode_range] = common_parts(x,...
                                    params,...
                                    init_shape)

    ode_step_size = params.L / (params.n_t - 1);
    ode_range     = 0 : ode_step_size : params.L;
    % Solve the BVP
    options = bvpset('Vectorized', 'on',...
                     'FJacobian' , @(s,p) odeJac(s,p, x, params),...
                     'BCJacobian', {diag([ ones(1,12) zeros(1,6)])
                                    diag([zeros(1,12)  ones(1,6)])}...
                     );
                 
    solinit = bvpinit(ode_range,...
                      @(r) init_shape(round(r/ode_step_size+1), :));

    sol = bvp4c(@( s, p)  ode_main(s,p, x, params),...
                @(xa,xb) BoundCond(xa,xb, params),...
                solinit,...
                options);

    % Compute dense output
    sxint = deval(sol, ode_range)';
    
    % xy distances (ignoring Z)
    dist = sqrt(sum(sxint(:,1:2).^2,2));
end





function sxint = final_shape(f_ly, params)

    sxint = common_parts(f_ly, params, params.initial_shape);
    
end



function dxds = ode_main(s,...
                         x,...
                         f_ly,...
                         params)

    % BOTTLENECK: spline() took up most computation time. But it was evaluated very often
    % for the same input [s] or equal spline representation [splines]. This means a lot of
    % work can be recycled between consecutive calls:
    persistent splines  f_y_prev  f_ly_prev  s_prev
% 
    % Recompute cubic splines, only when needed
    reeval = isempty(f_y_prev);
    if isempty(splines) || ~isequal(f_ly, f_ly_prev)
        splines = spline(0 : params.L/(numel(f_ly)-1) : params.L, f_ly);
        reeval  = true;
    end

    % Recompute value, only when needed
    f_y = f_y_prev;
    if reeval || ~isequal(s, s_prev)
        f_y = ppval_quick(splines, s); end

    f_y_prev  = f_y;
    s_prev    = s;
    f_ly_prev = f_ly;
    
    % look at syms_equations.m and symbolic_equations.m for jacobian and
    % vectorization
        
    C   = params.C;
    u_p = params.u_p;

    x4  = x( 4,:);   x5  = x( 5,:);    x6  = x( 6,:);     x7 = x( 7,:);
    x8  = x( 8,:);   x9  = x( 9,:);    x10 = x(10,:);    x11 = x(11,:);
    x12 = x(12,:);   x13 = x(13,:);    x14 = x(14,:);    x15 = x(15,:);
    x16 = x(16,:);   x17 = x(17,:);    x18 = x(18,:);

    f0  = params.rhoA  * params.g * (params.L - s);
    f1  = params.eg(1) * f0;
    f2  = params.eg(2) * f0;
    f3  = params.eg(3) * f0;

    f_yC = conj(f_y);

    G1 = params.EC_F(1) + conj(x16) + f1;      G4 = conj(u_p(3)) - conj(x15);
    G2 = params.EC_F(2) + conj(x17) + f2;      G5 = conj(u_p(2)) - conj(x14);
    G3 = params.EC_F(3) + conj(x18) + f3;      G6 = conj(u_p(1)) - conj(x13);

    dxds = [+x6
            +x9
            +x12
            +x5 .*x15 - x6 .*x14
            -x4 .*x15 + x6 .*x13
            +x4 .*x14 - x5 .*x13
            +x8 .*x15 - x9 .*x14
            -x7 .*x15 + x9 .*x13
            +x7 .*x14 - x8 .*x13
            +x11.*x15 - x12.*x14
            -x10.*x15 + x12.*x13
            +x10.*x14 - x11.*x13
            +(conj(x5).*G1 + conj(x8).*G2 + conj(x11).*G3 - C(2).*x15.*G5 + C(3).*x14.*G4)/C(1)
            -(conj(x4).*G1 + conj(x7).*G2 + conj(x10).*G3 - C(1).*x15.*G6 + C(3).*x13.*G4)/C(2)
            -(C(1)*x14.*G6 - C(2)*x13.*G5)/C(3)
            -x5 .*f_yC;
            -x8 .*f_yC;
            -x11.*f_yC];

end


function res = BoundCond(xa,...
                         xb,...
                         params)

    res = [ xa( 1:12) - [params.cyl_center; params.R(:)]
            xb(13:18) - [params.u_p; 0; 0; 0]
          ];
end



function dfdx = odeJac(s,...
                       x,...
                       f_ly,...
                       params)

    % BOTTLENECK: spline() took up most computation time. But it was evaluated very often
    % for the same input [s] or equal spline representation [splines]. This means a lot of
    % work can be recycled between consecutive calls:
%     reduced_basis = 0:params.L/(numel(f_ly)-1):params.L;
%     f_y = interp1(reduced_basis,f_ly,s);
    persistent splines  f_y_prev  f_ly_prev  s_prev

    reeval = isempty(f_y_prev);

    % Recompute cubic splines, only when needed
    if isempty(splines) || ~isequal(f_ly, f_ly_prev)
        splines = spline(0 : params.L/(numel(f_ly)-1) : params.L, f_ly);
        reeval  = true;
    end

    % Recompute value, only when needed
    f_y = f_y_prev;
    if reeval || ~isequal(s, s_prev)
        f_y = ppval_quick(splines, s); end

    f_y_prev  = f_y;
    s_prev    = s;
    f_ly_prev = f_ly;




    x4  = x( 4);       x5  = x( 5);    x6  = x( 6);     x7  = x( 7);
    x8  = x( 8);       x9  = x( 9);    x10 = x(10);     x11 = x(11);
    x12 = x(12);       x13 = x(13);    x14 = x(14);     x15 = x(15);
    x16 = x(16);       x17 = x(17);    x18 = x(18);

    c1 = params.C(1);
    c2 = params.C(2);
    c3 = params.C(3);

    f0 = params.g * params.rhoA * (params.L - s);
    f1 = params.eg(1)*f0;
    f2 = params.eg(2)*f0;
    f3 = params.eg(3)*f0;

    G1 = params.EC_F(1) + conj(x16) + f1;      G4 = c3 * (conj(params.u_p(3)) - conj(x15));
    G2 = params.EC_F(2) + conj(x17) + f2;      G5 = c1 * (conj(params.u_p(1)) - conj(x13));
    G3 = params.EC_F(3) + conj(x18) + f3;      G6 = c2 * (conj(params.u_p(2)) - conj(x14));

    f_yC = conj(f_y);


    dfdx = [0 0 0      0,       0,    1,       0,        0,    0,       0,       0,     0,                 0,                  0,                 0,            0,            0,             0
            0 0 0      0,       0,    0,       0,        0,    1,       0,       0,     0,                 0,                  0,                 0,            0,            0,             0
            0 0 0      0,       0,    0,       0,        0,    0,       0,       0,     1,                 0,                  0,                 0,            0,            0,             0
            0 0 0      0,     x15, -x14,       0,        0,    0,       0,       0,     0,                 0,                -x6,                x5,            0,            0,             0
            0 0 0   -x15,       0,  x13,       0,        0,    0,       0,       0,     0,                x6,                  0,               -x4,            0,            0,             0
            0 0 0    x14,    -x13,    0,       0,        0,    0,       0,       0,     0,               -x5,                 x4,                 0,            0,            0,             0
            0 0 0      0,       0,    0,       0,      x15, -x14,       0,       0,     0,                 0,                -x9,                x8,            0,            0,             0
            0 0 0      0,       0,    0,    -x15,        0,  x13,       0,       0,     0,                x9,                  0,               -x7,            0,            0,             0
            0 0 0      0,       0,    0,     x14,     -x13,    0,       0,       0,     0,               -x8,                 x7,                 0,            0,            0,             0
            0 0 0      0,       0,    0,       0,        0,    0,       0,     x15,  -x14,                 0,               -x12,               x11,            0,            0,             0
            0 0 0      0,       0,    0,       0,        0,    0,    -x15,       0,   x13,               x12,                  0,              -x10,            0,            0,             0
            0 0 0      0,       0,    0,       0,        0,    0,     x14,    -x13,     0,              -x11,                x10,                 0,            0,            0,             0
            0 0 0      0,  +G1/c1,    0,       0,   +G2/c1,    0,       0,  +G3/c1,     0,                 0,   (G4 + c2*x15)/c1, -(G6 + c3*x14)/c1,  conj(x5)/c1,  conj(x8)/c1,  conj(x11)/c1
            0 0 0 -G1/c2,       0,    0,  -G2/c2,        0,    0,  -G3/c2,       0,     0, -(G4 + c1*x15)/c2,                  0, +(G5 + c3*x13)/c2, -conj(x4)/c2, -conj(x7)/c2, -conj(x10)/c2
            0 0 0      0,       0,    0,       0,        0,    0,       0,       0,     0,  (G6 + c1*x14)/c3,  -(G5 + c2*x13)/c3,                 0,            0,            0,             0
            0 0 0      0,   -f_yC,    0,       0,        0,    0,       0,       0,     0,                 0,                  0,                 0,            0,            0,             0
            0 0 0      0,       0,    0,       0,    -f_yC,    0,       0,       0,     0,                 0,                  0,                 0,            0,            0,             0
            0 0 0      0,       0,    0,       0,        0,    0,       0,   -f_yC,     0,                 0,                  0,                 0,            0,            0,             0];

end


function X = ppval_quick(S, X)
% PPVAL   Compute the value of a cubic splines interpolant,
%         but then a bit faster than ppval() does it.
%
% See also spline.


    if isscalar(X) %(faster)
        % Find breakpoint
        br = find(X > S.breaks, 1, 'last');
        if isempty(br), br = 1; end

        % Compute spline interpolant
        X(:) = S.coefs(br,:) * (X - S.breaks(br)).^[3; 2; 1; 0];

    else % (more flexible)
        % Find breakpoint
        br = sum(bsxfun(@gt, X(:), S.breaks), 2);
        br(br==0) = 1;

        % Compute spline interpolant
        X(:) = sum(S.coefs(br,:) .* bsxfun(@power, X(:) - S.breaks(br).', [3 2 1 0]), 2);
    end

end % subfunction
