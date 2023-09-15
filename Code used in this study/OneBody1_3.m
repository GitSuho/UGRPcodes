%%% 1-body condition %%%
clear all
close all
clc
hold on;

%Initial consitions
% %circle
% x_init = [1, 0];
% v_init = [0, 3];
% M = 9;m = 1 ;G = 1;

%ellipse
x_init = [1.5, 0];
v_init = [0, 3];
M = 9;m = 1 ;G = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms('x_sym', [2, 2]);
assume(x_sym, 'real');
E = 1/2*m*norm(v_init)^2 - G*M*m/norm(x_init);
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
fprintf("Eccentricity : %d\n", l_const.^2*A_coeff/k_const);
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


dt = 1/100;
given_plotinterval = 0.1;

X_dt(1) = dt;
P_dt(1) = dt;

X(1, :) = init;
P(1, :) = init;

X_error(1, :) = [0, 0];
P_error(1, :) = [0, 0];

i = 1;
while (1)
    %fit each dt values
    X_dt(i+1) = Fit_dt(X_dt(i), given_plotinterval ,X(i,:), @geo);
    P_dt(i+1) = Fit_dt(P_dt(i), given_plotinterval ,P(i,:), @lag);
    
    %use Runge-Kutta method and predict next position
    X(i+1, :) = RK4(X(i, :), @geo, X_dt(i+1));
    P(i+1, :) = RK4(P(i, :), @lag, P_dt(i+1));
    
    if (Find_degree(X(end,1), X(end,2)) < pi)
        break;
    end

    %Error estimation according to analytic solution
    X_error(i+1, :) = Err_est(X(i+1,1), X(i+1,2));
    P_error(i+1, :) = Err_est(P(i+1,1), P(i+1,2));

    i = i+1;
end

%plot numerical solutions' trajectory
plot(X(:, 1), X(:, 2), 'bo');
plot(P(:, 1), P(:, 2), 'r^');
pause(1);

%write file
file_name = sprintf('norm_init[%2.2f,%2.2f,%2.2f,%2.2f]_dt%.5f_m[%2d,%2d].txt', init(1), init(2), init(3), init(4), dt, M, m);
wfile = fopen(file_name, 'w'); %file name은 initial value정하고 그 값에 맞춰 하기. 파일의 첫째줄은 그대로

format1 = 'Error of %s | mean = %f, max = %f, min = %f\n';
format2 = '%3d||%f | %f | %f || %f | %f | %f \n';
format3 = '   ||      arc domain geodesic      ||           newtonian           \nnum||degree   | error    | dt       || degree   | error    | dt       \n'; 
fprintf(wfile, "%s\nEccentricity : %f \n\n", file_name , l_const.^2*A_coeff/k_const);
fprintf(wfile, format1, 'arcl_geod', mean(X_error(2:length(X_error), 2)), max(X_error(2:length(X_error), 2)), min(X_error(2:length(X_error),2 )) );
fprintf(wfile, format1, 'newtonian', mean(P_error(2:length(P_error), 2)), max(P_error(2:length(P_error), 2)), min(P_error(2:length(P_error),2 )) );
fprintf(wfile, '\n');

fprintf(wfile, format3);
for i = 1:length(X_error)
    fprintf(wfile, format2, i ,X_error(i,:) ,X_dt(i) ,  P_error(i,:), P_dt(i));
end


fclose(wfile);
hold off;
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
%newtonian equation
function slope = lag(x)
    global EL;
    slope = zeros(size(x));
    for i = 1:2
        slope(i) = x(i+2);
        slope(i+2) = EL{i}(x);
    end
end
%Euler 1st Order
function next_pos = EU1(curr_pos, func, dt)
    next_pos = curr_pos + func(curr_pos)*dt;
end
%Heun 2nd Order
function next_pos = HE2(curr_pos, func, dt)
    k1 = func(curr_pos);
    k2 = func(curr_pos + k1*dt);
    next_pos = curr_pos + (k1+k2)*dt/2;
end
%Kutta 3rd Order
function next_pos = KU3(curr_pos, func, dt)
    k1 = func(curr_pos);
    k2 = func(curr_pos + k1*dt/2);
    k3 = func(curr_pos - k1*dt + 2*k2*dt);
    next_pos = curr_pos + (k1 + 4*k2 + k3)*dt/6;
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
%find a dt value that make given plot interval
function result_dt = Fit_dt(dt, fit_val, curr_pos, func)
    accuracy = 0.00000001;
    next_pos = RK4(curr_pos,func, dt);
    disp  = sqrt((next_pos(1) - curr_pos(1)).^2 + (next_pos(2) - curr_pos(2)).^2);
    while( abs(disp - fit_val) > accuracy)
        dt = fit_val*dt/disp;
        next_pos = RK4(curr_pos,func, dt);
        disp  = sqrt((next_pos(1) - curr_pos(1)).^2 + (next_pos(2) - curr_pos(2)).^2);
    end
    result_dt = dt;
end
%calculate error when x-y coordinate is given
function result_err = Err_est(x, y)
    theta = Find_degree(x, y);
    simple_err = abs(OrbEqu_r(theta) - norm(x, y));
    result_err = [theta simple_err];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%