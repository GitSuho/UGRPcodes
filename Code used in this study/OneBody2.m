%%% 1-body condition %%%
clear all;

E_list = [-0.05, -0.1, -0.15, -0.20, -0.25];
ecc_list = [0.9, 0.75, 0.6 , 0.45 , 0.3 , 0.15];
numerical_method = ["EU1","HE2","KU3","RK4"];

for ord = 1:1
    for ii = 1:1
        for jj = 1:length(ecc_list)
close all;
hold on;

E = E_list(ii);
ecc = ecc_list(jj);
nume_name = numerical_method(ord);

filename = sprintf("OneBody2_E%1.2f_ecc%1.1f_%s.txt", E, ecc, nume_name);
figurename = sprintf("OneBody2_E%1.2f_ecc%1.1f_%s.pdf", E, ecc, nume_name);
fig1 = figure(1);

fprintf(filename);

M = 1;m = 1 ;G = 1;
v_init = [0, sqrt(-E/(-m/2+m/(ecc+1)))];
x_init = [(ecc+1)*G*M/(v_init(2).^2), 0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms('x_sym', [2, 2]);
assume(x_sym, 'real');
E = 1/2*m*norm(v_init)^2 - G*M*m/norm(x_init);
T = E + G*M*m/norm(x_sym(1,:));


init = [x_init v_init];
init = init(:);
%Geodesic solution which domain is arclengthP
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
xlabel("x"); ylabel("y");

dt = 1/50;

X(1, :) = init;
P(1, :) = init;

X_error(1, :) = [0, 0, 0];
P_error(1, :) = [0, 0, 0];


%write file
wfile = fopen(filename, 'w');

fprintf(wfile, filename);
fprintf(wfile, ['\ninitial value : [x, y, v_x, v_y] = [%2.2f,%2.2f,%2.2f,%2.2f] , ' ...
    'numerical interval dt = %.5f , [M, m] = [%2d,%2d]\n\n'], init(1), init(2), init(3), init(4), dt, M, m);
fprintf(wfile, "num|geodesic|\nnum|newtonian\nnum|degree, x, y, simple error, relative error\n");

format = "%3d|%f,%f,%f,%f,%f\n";

F_count = 0;
fprintf(wfile, format, 1, X_error(1,1), X(1,1), X(1,2), X_error(1,2), X_error(1,3));
while(1)
    for i = 2:101
        %use Runge-Kutta method and predict next position
        X(i, :) = RKn([OrbEqu(X_error(i-1,1)), X(i-1,3:4)], @geo, dt, ord);    
        %Error estimation according to analytic solution
        X_error(i, :) = Err_est(X(i,1), X(i,2));    

        if (X_error(i, 1) < X_error(i-1, 1))
            break;
        end
    end
    
    plot(X(1:i,1), X(1:i,2), 'r-o');
    for j = 2:i
    fprintf(wfile, format, j + F_count*100, X_error(j,1), X(j,1), X(j,2), X_error(j,2), X_error(j,3));
    end
    
    if(X_error(101, 1) < X_error(100, 1) || i < 101)
        break;
    end

    X(1, :) = X(101, :);
    X_error(1, :) = X_error(101, :);  
    F_count = F_count+1;
end


F_count = 0;
fprintf(wfile, format, 1, P_error(1,1), P(1,1), P(1,2), P_error(1,2), P_error(1,3));
while(1)
    for i = 2:101
        %use Runge-Kutta method and predict next position
        P(i, :) = RKn([OrbEqu(P_error(i-1,1)), P(i-1,3:4)], @lag, dt, ord);    
        %Error estimation according to analytic solution
        P_error(i, :) = Err_est(P(i,1), P(i,2));    

        if (P_error(i, 1) < P_error(i-1, 1))
            break;
        end
    end
    
    plot(P(1:i,1), P(1:i,2), 'b-+');
    for j = 2:i
    fprintf(wfile, format, j + F_count*100, P_error(j,1), P(j,1), P(j,2), P_error(j,2), P_error(j,3));
    end
    
    if(P_error(101, 1) < P_error(100, 1) || i < 101)
        break;
    end

    P(1, :) = P(101, :);
    P_error(1, :) = P_error(101, :);  
    F_count = F_count+1;
end


fclose(wfile);
hold off;
exportgraphics(fig1, figurename);


fprintf('program end\n');

        end
    end
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
        case 2
            k1 = func(curr_pos);
            k2 = func(curr_pos + k1*dt);
            next_pos = curr_pos + (k1+k2)*dt/2;
        case 1
            next_pos = curr_pos + func(curr_pos)*dt;
    end
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
function result_err = Err_est(x, y)
    theta = Find_degree(x, y);
    simple_err = abs(OrbEqu_r(theta) - norm(x, y));
    relative_err = simple_err / OrbEqu_r(theta) * 100;
    result_err = [theta, simple_err, relative_err];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%