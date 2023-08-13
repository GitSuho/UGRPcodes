%1-body by Cho


syms('x_sym', [2, 2]);
%[x1_1 x1_2  = [x1 x2
% x2_1 x2_2]    v1 v2]
assume(x_sym, 'real');

x_init = [1, 2];
v_init = [-0.5, -2];
M = 1;m = 10;
G = 1;


E = 1/2*m*norm(v_init)^2 - G*M*m/norm(x_init);
T = E + G*M*m/norm(x_sym(1,:));
global F;
F{2} = [];
f{2} = [];
for i = 1:2
    f{i} = -(diff(T, x_sym(1,1))*x_sym(2,1) + diff(T, x_sym(1,2))*x_sym(2,2))*x_sym(2,i) ...
       + 1/(2*T)*diff(T, x_sym(1,i))*norm(x_sym(2,:))^2;
    F{i} = @(x) eval(subs(f{i}, x_sym, reshape(x,[2,2])));
end
global EL;
EL{2} = [];
el{2} = [];
for i = 1:2
    el{i} = -G*M*x_sym(1, i)/(norm(x_sym(1,:).^1.5));
    EL{i} = @(x) eval(subs(el{i}, x_sym, reshape(x,[2,2])));
end


init = [x_init v_init];
init = init(:);


dt = 1/100;
cycle = 10;
clear("X")
clear("P")
X(1, :) = init;
P(1, :) = init;

plot(0, 0, 'blacko');

while(1)
    for i = 1:cycle
        % X(i+1, :) = RK4(X(i, :), @g, dt);
        P(i+1, :) = RK4(P(i, :), @n, dt);
    end

    hold off;
    hold on;

    % plot(X(:, 1), X(:, 2), 'b-');
    plot(P(:, 1), P(:, 2), 'r:');

    % X(1, :) = X(cycle+1, :);
    P(1, :) = P(cycle+1, :);
    pause(0.1);
end


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






