%%% 1-body condition %%%
clear all
close all
clc

%Initial consitions


E_const = -0.5;
ecc_lis = [0 0.25 0.5 0.75];

for ind = 1:4
    ecc = ecc_lis(ind);
    close all
    hold on;

M = 1;m = 1 ;G = 1;
v_init = [0, sqrt(E_const/(-m/2+m/(ecc+1)))];
x_init = [(ecc+1)*G*M/(v_init(2).^2), 0];

% fprintf("initial v : %f, x : %f\n", v_init(2), x_init(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms('x_sym', [2, 2]);
assume(x_sym, 'real');
E = 1/2*m*norm(v_init)^2 - G*M*m/norm(x_init);

fprintf("E : %f\n", E);

T = E + G*M*m/norm(x_sym(1,:));
init = [x_init v_init];
init = init(:);
%Geodesic Solution which domain is arclength
global GE;
GE{2} = [];
ge{2} = [];
for i = 1:2
    ge{i} = -(diff(T, x_sym(1,1))*x_sym(2,1) + diff(T, x_sym(1,2))*x_sym(2,2))*x_sym(2,i) ...
       + 1/(2*T)*diff(T, x_sym(1,i))*norm(x_sym(2,:))^2;
    GE{i} = @(x) eval(subs(ge{i}, x_sym, reshape(x,[2,2])'));
end
%Euler-Lagrange solution
global EL;
EL{2} = [];
el{2} = [];
for i = 1:2
    el{i} = -G*M*x_sym(1, i)/(norm(x_sym(1,:))^3);
    EL{i} = @(x) eval(subs(el{i}, x_sym, reshape(x,[2,2])'));
end

%Orbital equation's constants
global l_const; global A_coeff; global k_const;
k_const = G*M;
l_const = norm(x_init)*norm(v_init)*sin(acos((norm(x_init).^2+norm(v_init).^2-norm(x_init-v_init).^2) ...
    /(2*norm(x_init)*norm(v_init))));%r*v*sin(theta_0)
A_coeff =  (1/norm(x_init) - (k_const/(l_const.^2)))/(x_init(1)/norm(x_init));

% fprintf("Eccentricity : %f\n", l_const.^2*A_coeff/k_const);
% fprintf("Eccentricity : %f\n", (((x_init(1)*v_init(2)-x_init(2)*v_init(1)).^2/(G*M*norm(x_init)))-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot an analytic Solution's trajectory
count = 10000;
orbit_coordinate(1, :) = [0 0];
theta_arr = linspace(0, 2*pi, count);
for i = 1:count
    orbit_coordinate(i, :) = OrbEqu(theta_arr(i));
end
plot(orbit_coordinate(:, 1), orbit_coordinate(:, 2), 'black-');
plot(0, 0, 'blacko');
pause(0.1);


dt = 1/2;

X(1, :) = init;
P(1, :) = init;

X_error(1, :) = [0, 0];
P_error(1, :) = [0, 0];

file_name = sprintf('Kims_request_e%2.2f_E%2.2f_init[%2.2f,%2.2f,%2.2f,%2.2f]_dt%.5f_m[%2d,%2d].txt', l_const.^2*A_coeff/k_const, E, init(1), init(2), init(3), init(4), dt, M, m);
wfile = fopen(file_name, 'w');

for period = 0:0
    fprintf(wfile, "cycle %d\n", period+1);
    i = 1;
    while(1)
        plot(X(i,1), X(i,2), 'ro');
        pause(0.1);


        X(i+1, :) = RK4(X(i, :), @geo, dt);
        if(Find_degree(X(i+1,1), X(i+1,2)) - Find_degree(X(i,1), X(i,2)) < 0)
            break;
        end
        X_error(i+1, :) = Err_est(X(i+1,1), X(i+1,2), period);
        i = i+1;
    end
    i = 1;
    while(1)
        plot(P(i,1), P(i,2), 'b+');
        pause(0.1);

        P(i+1, :) = RK4(P(i, :), @lag, dt);
        if(Find_degree(P(i+1,1), P(i+1,2)) - Find_degree(P(i,1), P(i,2)) < 0)
            break;
        end
        P_error(i+1, :) = Err_est(P(i+1,1), P(i+1,2), period);
        i = i +1;
    end
    var1 = X(end,:);
    var2 = P(end,:);
    
    X(end,:) = [];
    P(end,:) = [];
    

    if (length(X) > length(P))
        len = length(P);
        for j = (1:length(X)-length(P))
            P_error(len+j, :) = [0 0];
            P(len+j, :) = [0 0 0 0];
        end
    elseif (length(X) < length(P))
        len = length(X);
        for j = 1:(length(P)-length(X))
            X_error(len+j, :) = [0 0];
            X(len+j, :) = [0 0 0 0];
        end
    end

    for i = 1:length(X)
        fprintf(wfile, "%f,%f,%f,%f|%f,%f,%f,%f\n", X(i,1), X(i,2), X_error(i,1), X_error(i,2), P(i,1), P(i,2), P_error(i,1), P_error(i,2));
        %%% x1, x2, x_theta, x_error | p1, p2, p_theta, p_error
    end

    clear X_error P_error X P
    X_error(1,:) = Err_est(var1(1), var1(2), period+1);
    P_error(1,:) = Err_est(var2(1), var2(2), period+1);
    X(1,:) = var1;
    P(1,:) = var2;
end


fclose(wfile);
% pause(3);
% hold off;
fprintf('program end %d\n', ind);

end

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
%newtonian equation
function slope = lag(x)
    global EL;
    slope = zeros(size(x));
    for i = 1:2
        slope(i) = x(i+2);
        slope(i+2) = EL{i}(x);
    end
end
%Runge-Kutta 4th Order
function next_pos = RK4(curr_pos, func, dt)
    k1 = func(curr_pos);
    k2 = func(curr_pos + k1*dt/2); 
    k3 = func(curr_pos + k2*dt/2); 
    k4 = func(curr_pos + k3*dt); 
    next_pos = curr_pos + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
end
%calculate point of a moving particle when theta is given
function xy_coor = OrbEqu(theta)
    global l_const; global A_coeff; global k_const;
    r = 1/(A_coeff*cos(theta)+(k_const/(l_const.^2)));
    xy_coor = [r*cos(theta), r*sin(theta)];
end
%calculate the r which is interval of two particles
function r = OrbEqu_r(theta)
    global l_const; global A_coeff; global k_const;
    r = 1/(A_coeff*cos(theta)+(k_const/(l_const.^2)));
end
%calculate a degree when x-y coordinate given
function degree = Find_degree(x, y)
    if (y >= 0 )
        degree = acos( x/sqrt(x.^2 + y.^2));
    else
        degree = 2*pi - acos( x/sqrt(x.^2 + y.^2));
    end
end
%calculate error when x-y coordinate is given
function result_err = Err_est(x, y, period)
    theta = Find_degree(x, y);
    simple_err = abs(OrbEqu_r(theta) - norm(x, y));
    result_err = [theta+2*pi*period simple_err];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%