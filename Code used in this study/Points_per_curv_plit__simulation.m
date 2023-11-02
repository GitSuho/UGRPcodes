%%% 1-body condition %%%
clear all;

curv_lis = [0.001:0.05:2.001];
plit_lis = [0.001:0.05:2.001];

ord = 2;
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for curv = curv_lis
    % for plot_interval = plit_lis

        curv = 1;
        curv

        M = 0.1;m = 1 ;G = 1;
        v_init = [0, 1];
        x_init = [1/curv, 0];

        % ecc = 0;
        % E = -0.1;
        % 
        % M = 9;m = 1 ;G = 1;
        % v_init = [0, sqrt(-E/(-m/2+m/(ecc+1)))];
        % x_init = [(ecc+1)*G*M/(v_init(2).^2), 0];
        % curv=1/x_init(1);
        % curv

        syms('x_sym', [2, 2]);
        assume(x_sym, 'real');
        E = 1/2*m*norm(v_init)^2 - G*M*m/norm(x_init);
        T = E + G*M*m/norm(x_sym(1,:));
        init = [x_init v_init];
        init = init(:);
        %Geodesic solution which domain is arclength
        global GE;
        GE{2} = [];
        ge{2} = [];
        for i = 1:2
            ge{i} = -(diff(T, x_sym(1,1))*x_sym(2,1) + diff(T, x_sym(1,2))*x_sym(2,2))*x_sym(2,i) ...
               + 1/(2*T)*diff(T, x_sym(1,i))*norm(x_sym(2,:))^2;
            GE{i} = @(x) eval(subs(ge{i}, x_sym, reshape(x,[2,2])'));
        end
        %Orbital equation's constants
        global l_const; global A_coeff; global k_const;
        k_const = G*M;
        l_const = norm(x_init)*norm(v_init)*sin(acos((norm(x_init).^2+norm(v_init).^2-norm(x_init-v_init).^2) ...
            /(2*norm(x_init)*norm(v_init))));%r*v*sin(theta_0)
        A_coeff =  (1/norm(x_init) - (k_const/(l_const.^2)))/(x_init(1)/norm(x_init));
        
        %draw exact trajectory
        count = 10000;
        orbit_coordinate(1, :) = [0 0];
        theta_arr = linspace(0, 2*pi, count);
        for i = 1:count
            orbit_coordinate(i, :) = OrbEqu(theta_arr(i));
        end
        plot(orbit_coordinate(:, 1), orbit_coordinate(:, 2), 'black-');
        plot(0, 0, 'blacko'); 
        xlabel("x"); ylabel("y");

        
        plot_interval = 0.5;

        dt = 1/50;       
        dt = Fit_dt(dt, plot_interval ,init, @geo, ord);
        X = RKn(init, @geo, dt, ord);
        plot([init(1), X(1)], [init(2), X(2)], 'r:^' );


        % plot_interval = norm(init)*dt;
        X_theo = TheoPredic_n(init, plot_interval, ord);
        plot([init(1), X_theo(1)], [init(2), X_theo(2)], 'b:o' );
        


        
        % X_error = Err_est(X(1), X(2), 1/curv);
    % end
% end
fprintf('program end\n');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%geodesic equation
function dxdt = geo(x)
    global GE;
    dxdt = zeros(size(x));
    for i = 1:2
        dxdt(i) = x(i+2);
        dxdt(i+2) = GE{i}(x);
    end
end
%calculate point of a moving particle when theta is given
function xy_coor = OrbEqu(theta)
    global l_const; global A_coeff; global k_const;
    r = 1/(A_coeff*cos(theta)+(k_const/(l_const.^2)));
    xy_coor = [r*cos(theta), r*sin(theta)];
end
function xyv1v2 = OrbEqu_with_velocity(theta)
    global l_const; global A_coeff; global k_const;
    r = 1/(A_coeff*cos(theta)+(k_const/(l_const.^2)));
    r_p = r.^2 * A_coeff*sin(theta);%derivative r by theta
    % assume that [ d{\theta}/dt = 1 ]
    xyv1v2 = [r*cos(theta), r*sin(theta) , r_p*cos(theta)-r*sin(theta) , r_p*sin(theta)+r*cos(theta)];
end
%Runge-Kutta nth Order
function next_pos = RKn(curr_pos, func, dt, order)
    switch order
        case 4
            k1 = func(curr_pos);
            k2 = func(curr_pos + k1*dt/2); 
            k3 = func(curr_pos + k2*dt/2); 
            k4 = func(curr_pos + k3*dt); 
            next_pos = curr_pos + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
        case 3
            k1 = func(curr_pos);
            k2 = func(curr_pos + k1*dt/2);
            k3 = func(curr_pos - k1*dt + 2*k2*dt);
            next_pos = curr_pos + (k1 + 4*k2 + k3)*dt/6;
        case 2      %임시적으로 Euler mid point method로 바꿨음.
            k1 = func(curr_pos);
            k2 = func(curr_pos + k1*dt/2);
            next_pos = curr_pos + (k2)*dt;
        case 1
            next_pos = curr_pos + func(curr_pos)*dt;
    end
end
%Theoritical predict nth Order
function next_pos = TheoPredic_n(curr_pos, d, order)
    switch order
        case 4

        case 3

        case 2      
            r = curr_pos(1);
            next_pos = [r - d.^2./sqrt(4.*r.^2 + d.^2) , 2.*r.*d./sqrt(4.*r.^2 + d.^2)];
        case 1

    end
end
%calculate the r which is interval of two particles
function r = OrbEqu_r(theta)
    global l_const; global A_coeff; global k_const;
    r = 1/(A_coeff*cos(theta)+(k_const/(l_const.^2)));
end
%calculate a degree when x-y coordinate given
function degree = Find_degree(x, y)
    degree = mod (atan2(y, x), 2*pi);
end
%find a curvature at given angle
function curvature = Find_Curv(angle_0, ecc)
    global a_const; global b_const;
    angle = 2*atan(sqrt((1-ecc)/(1+ecc))*tan(angle_0/2));
    curvature = a_const*b_const/(((a_const*sin(angle)).^2 + (b_const*cos(angle)).^2).^(1.5));
end
%find a dt value that make given plot interval
function result_dt = Fit_dt(dt, fit_val, curr_pos, func, ord)
    accuracy = 0.00000001;
    next_pos = RKn(curr_pos,func, dt, ord);
    disp  = sqrt((next_pos(1) - curr_pos(1)).^2 + (next_pos(2) - curr_pos(2)).^2);
    while( abs(disp - fit_val) > accuracy)
        dt = fit_val*dt/disp;
        next_pos = RKn(curr_pos,func, dt, ord);
        disp  = sqrt((next_pos(1) - curr_pos(1)).^2 + (next_pos(2) - curr_pos(2)).^2);
    end
    result_dt = dt;
end
%calculate error when x-y coordinate is given
function result_err = Err_est(x, y, r_0)
    % rr = OrbEqu_r(Find_degree(x, y));
    result_err = 100 * abs(r_0 - sqrt(x.^2 + y.^2)) / r_0; %relative_err
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%