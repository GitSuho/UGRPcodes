%%% 1-body condition %%%
close all
clc

%Initial consitions
E_const = 0.5;
ecc = 0.5;


hold on;

M = 1;m = 1 ;G = 1;
v_init = [0, sqrt(E_const/(-m/2+m/(ecc+1)))];
x_init = [(ecc+1)*G*M/(v_init(2).^2), 0];

fprintf("initial v : %f, x : %f\n", v_init(2), x_init(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms('x_sym', [2, 2]);
assume(x_sym, 'real');
E = 1/2*m*norm(v_init)^2 - G*M*m/norm(x_init);

fprintf("E : %f\n", E);

%Orbital equation's constants
global l_const; global A_coeff; global k_const;
k_const = G*M;
l_const = norm(x_init)*norm(v_init)*sin(acos((norm(x_init).^2+norm(v_init).^2-norm(x_init-v_init).^2) ...
    /(2*norm(x_init)*norm(v_init))));%r*v*sin(theta_0)
A_coeff =  (1/norm(x_init) - (k_const/(l_const.^2)))/(x_init(1)/norm(x_init));
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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