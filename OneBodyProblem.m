%%% 1-body condition %%%

%Initial consitions
x_init = [1, 0];
v_init = [0, 3];
M = 9;m = 1 ;G = 1;

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


%Geodesic Solution which domain is time
global GE_t;
GE_t{2} = [];
ge_t{2} = [];
for i = 1:2
    ge_t{i} = (-G*M*m*dot(x_sym(1,:),x_sym(2,:))/(norm(x_sym(1,1:2))^3))/T*x_sym(2, i) ...
    -(1/T)*(diff(T, x_sym(1, 1))*x_sym(2, 1) + diff(T, x_sym(1, 2))*x_sym(2, 2))*x_sym(2, i) ... 
    + 1/(2*T)*diff(T, x_sym(1, i))*(x_sym(2, 1)^2 + x_sym(2, 2)^2);
    GE_t{i} = @(x) eval(subs(ge_t{i}, x_sym, reshape(x,[2,2])'));
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

hold on;

%Analytic Solution
count = 10000;
orbit_coordinate(1, :) = [0 0];
theta_arr = linspace(0, 2*pi, count);
for i = 1:count
    orbit_coordinate(i, :) = OrbEqu(theta_arr(i));
end
plot(orbit_coordinate(:, 1), orbit_coordinate(:, 2), 'black-');
plot(0, 0, 'blacko');


dt = 1/100;
cycle = 10;
clear("X")
clear("X_t")
clear("P")
X(1, :) = init;
X_t(1, :) = init;
P(1, :) = init;

% %plot numerical solutions
% while(1)
%     for i = 1:cycle
%         X(i+1, :) = RK4(X(i, :), @geo, dt);
%         X_t(i+1, :) = RK4(X_t(i, :), @geo_t, dt);
%         P(i+1, :) = RK4(P(i, :), @lag, dt);
%     end
%     plot(X(:, 1), X(:, 2), 'b--o');
%     plot(X_t(:, 1), X_t(:, 2), 'g--p');
%     plot(P(:, 1), P(:, 2), 'r--^');
% 
%     X(1, :) = X(cycle+1, :);
%     X_t(1, :) = X_t(cycle+1, :);
%     P(1, :) = P(cycle+1, :);
%     pause(0.1);
% end

%calculate error values
X_count = 104*100;
X_t_count = 104*100;
P_count = 104*100;
X_dt = 1/50; 
X_t_dt = 1/50; 
P_dt = 1/50;

clear("X_error")
clear("X_t_error")
clear("P_error")
X_error(1, :) = [0, 0];
X_t_error(1, :) = [0, 0];
P_error(1, :) = [0, 0];

for i = 1:X_count-1
    X(i+1, :) = RK4(X(i, :), @geo, X_dt);
    % foo = norm(X(i+1, 3:4));
    % X(i+1, 3:4) = 20*X(i+1, 3:4)/foo;
end
for i = 1:X_t_count-1
    X_t(i+1, :) = RK4(X_t(i, :), @geo_t, X_t_dt);
end
for i = 1:P_count-1
    P(i+1, :) = RK4(P(i, :), @lag, P_dt);
    % foo = norm(P(i+1, 3:4));
    % P(i+1, 3:4) = P(i+1, 3:4)/foo;
end

plot(X(:, 1), X(:, 2), 'ro');
plot(X_t(:, 1), X_t(:, 2), 'gp');
plot(P(:, 1), P(:, 2), 'b^');

for i = 1:X_count
    foo_theta = Find_degree(X(i,1), X(i,2));
    var_error = abs(OrbEqu_r(foo_theta) - norm(X(i, 1:2)));
    X_error(i, :) = [foo_theta var_error];
end
for i = 1:X_t_count
    foo_theta = Find_degree(X_t(i,1), X_t(i,2));
    var_error = abs(OrbEqu_r(foo_theta) - norm(X_t(i, 1:2)));
    X_t_error(i, :) = [foo_theta var_error];
end
for i = 1:P_count
    foo_theta = Find_degree(P(i, 1), P(i, 2));
    var_error = abs(OrbEqu_r(foo_theta) - norm(P(i, 1:2)));
    P_error(i, :) = [foo_theta var_error];
end

fprintf("arc domain geo_error - mean : %d , max : %d , min : %d\n" , mean(X_error(2:X_count, 2)), max(X_error(2:X_count, 2)), min(X_error(2:X_count,2 ))); 
fprintf("time domain geo_error - mean : %d , max : %d , min : %d\n" , mean(X_t_error(2:X_t_count, 2)), max(X_t_error(2:X_t_count, 2)), min(X_t_error(2:X_t_count,2 ))); 
fprintf("newtonian_error - mean : %d , max : %d , min : %d\n" , mean(P_error(2:P_count, 2)), max(P_error(2:P_count, 2)), min(P_error(2:P_count,2 ))); 
xlswrite('arc domain geo_error', X_error);
xlswrite('time domain geo_error', X_t_error);
xlswrite('newtonian_error', P_error);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdt = geo(x)
    global GE;
    dxdt = zeros(size(x));
    for i = 1:2
        dxdt(i) = x(i+2);
        dxdt(i+2) = GE{i}(x);
        % dxdt(i) = 0.1*x(i+2)/norm(x(3:4));
        % dxdt(i+2) = 0.1*GE{i}(x);
    end
end

function dxdt = geo_t(x)
    global GE_t;
    dxdt = zeros(size(x));
    for i = 1:2
        dxdt(i) = x(i+2);
        dxdt(i+2) = GE_t{i}(x);
        % dxdt(i) = 0.1*x(i+2)/norm(x(3:4));
        % dxdt(i+2) = 0.1*GE_t{i}(x);
    end
end

function slope = lag(x)
    global EL;
    slope = zeros(size(x));
    for i = 1:2
        slope(i) = x(i+2);
        slope(i+2) = EL{i}(x);
        % slope(i) = 0.1*x(i+2)/norm(x(3:4));
        % slope(i+2) = 0.1*EL{i}(x);
    end
end

function next_pos = RK4(curr_pos, func, dt)
    k1 = func(curr_pos);
    k2 = func(curr_pos + k1*dt/2); 
    k3 = func(curr_pos + k2*dt/2); 
    k4 = func(curr_pos + k3*dt); 
    next_pos = curr_pos + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
end

function xy_coor = OrbEqu(theta)
    global l_const; global A_coeff; global k_const;
    r = 1/(A_coeff*cos(theta)+(k_const/(l_const.^2)));
    xy_coor = [r*cos(theta), r*sin(theta)];
end

function r = OrbEqu_r(theta)
    global l_const; global A_coeff; global k_const;
    r = 1/(A_coeff*cos(theta)+(k_const/(l_const.^2)));
end

function degree = Find_degree(x, y)
    if (y >= 0 )
        degree = acos( x/sqrt(x.^2 + y.^2));
    else
        degree = 2*pi - acos( x/sqrt(x.^2 + y.^2));
    end
end