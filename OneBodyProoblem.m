%1-body condition. 


syms('x_sym', [2, 2]);
assume(x_sym, 'real');

%hyperbola?
% x_init = [1, 0];
% v_init = [0, 3.5];
% M = 5;m = 2;G = 1;

x_init = [5, 1];
v_init = [1, 2];
M = 10;m = 2;G = 1;


E = 1/2*m*norm(v_init)^2 - G*M*m/norm(x_init);
T = E + G*M*m/norm(x_sym(1,:));

%Geodesic Solution which domain is arclength
global F;
F{2} = [];
f{2} = [];
for i = 1:2
    f{i} = -(diff(T, x_sym(1,1))*x_sym(2,1) + diff(T, x_sym(1,2))*x_sym(2,2))*x_sym(2,i) ...
       + 1/(2*T)*diff(T, x_sym(1,i))*(x_sym(2,1).^2+x_sym(2,2).^2);
    F{i} = @(x) eval(subs(f{i}, x_sym, reshape(x,[2,2])));
end

%Euler-Lagrange solution
global EL;
EL{2} = [];
el{2} = [];
for i = 1:2
    el{i} = -G*M*x_sym(1, i)/(norm(x_sym(1,:)).^(3/2));
    EL{i} = @(x) eval(subs(el{i}, x_sym, reshape(x,[2,2])));
end

%orbital equation's constants
a = norm(x_init);
b = norm(v_init);
c = norm(x_init-v_init);
global l_const; global A_coeff; global k_const;
k_const = G*M;
l_const = norm(x_init)*norm(v_init)*sin(acos((a.^2+b.^2-c.^2)/(2*a*b)));%r*v*sin(theta_0)
A_coeff =  (1/norm(x_init) - (k_const/(l_const.^2)))/(x_init(1)/norm(x_init));


init = [x_init v_init];
init = init(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% count = 10000;
% orbit_coordinate(1, :) = [0 0];
% theta_arr = linspace(0, 2*pi, count);
% for i = 1:count
%     orbit_coordinate(i, :) = OrbEqu(theta_arr(i));
% end
% plot(orbit_coordinate(:, 1), orbit_coordinate(:, 2), 'black-');


% F1 = -(diff(T, x_sym(1,1))*x_sym(2,1) + diff(T, x_sym(1,2))*x_sym(2,2))*x_sym(2,1) ...
%        + 1/(2*T)*diff(T, x_sym(1,1))*norm(x_sym(2,:))^2;
% F2 = -(diff(T, x_sym(1,1))*x_sym(2,1) + diff(T, x_sym(1,2))*x_sym(2,2))*x_sym(2,2) ...
%        + 1/(2*T)*diff(T, x_sym(1,2))*norm(x_sym(2,:))^2;
% f1 = @(t) eval(subs(f{1}, x_sym, reshape(t,[2,2])'));
% f2 = @(t) eval(subs(f{2}, x_sym, reshape(t,[2,2])'));
% % [x1 x2 v1 v2] -> [x1 x2
% %                   v1 v2]
% g = @(t,x) [x(3);
%             x(4);
%             f1(x);
%             f2(x)];
% 
% [~, x] = ode45(g, linspace(0, 1, 100), [x_init v_init]);
% plot(x(:,1),x(:,2),'linewidth',3);



dt = 1/100;
cycle = 10;
clear("X")
clear("P")
X(1, :) = init;
P(1, :) = init;

plot(0, 0, 'blacko');

while(1)
    for i = 1:cycle
        X(i+1, :) = RK4(X(i, :), @g, dt);
        % P(i+1, :) = RK4(P(i, :), @n, dt);
    end

    hold off;
    hold on;

    plot(X(:, 1), X(:, 2), 'b-');
    % plot(P(:, 1), P(:, 2), 'r:');

    X(1, :) = X(cycle+1, :);
    % P(1, :) = P(cycle+1, :);
    pause(0.1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdt = g(x)
    global F;
    dxdt = zeros(size(x));
    for i = 1:2
        dxdt(i) = x(i+2);
        dxdt(i+2) = F{i}(x);
    end
end

function slope = n(x)
    global EL;
    slope = zeros(size(x));
    for i = 1:2
        slope(i) = x(i+2);
        slope(i+2) = EL{i}(x);
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





