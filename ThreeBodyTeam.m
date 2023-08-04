x_init = [-1/2 -1/3; 1/2 -1/3; 0 2/3]; 
v_init = [1/3 -2/3; -1/3 2/3; -1 0]; 
m = [1 1 1];
G = 9.8;
init = [x_init v_init]';
init = init(:)';


E = 0;
for i = 0:2
    j = mod(i,3)+1;
    k = mod(i+1,3)+1;
    E = E + 1/2*m(j)*norm(v_init(j,:))^2 ...
        - G*m(j)*m(k)/norm(x_init(j,:)-x_init(k,:));
end

syms('x_sym', [3 4]);
assume(x_sym, 'real')

% T = E - V
T_sym = E;
for i = 0:2
    j = mod(i,3)+1;
    k = mod(i+1,3)+1;
    T_sym = T_sym + G * m(j) * m(k) / norm(x_sym(j,1:2)-x_sym(k,1:2));
end

T_sym_diff = 0;
for i = 0:2
    j = mod(i,3)+1;
    k = mod(i+1,3)+1;
    T_sym_diff = T_sym_diff + G * m(j) * m(k) * dot(x_sym(j,1:2)-x_sym(k,1:2),x_sym(j,3:4)-x_sym(k,3:4)) ...
                                / norm(x_sym(j,1:2)-x_sym(k,1:2))^2 ;
end

global F;
F{3,2} = [];
f{3,2} = [];
for i = 1:3
    for j = 1:2
         f{i,j} = T_sym_diff / T_sym * x_sym(i,j+2);
         for k = 1:3
             for l = 1:2
                f{i,j} = f{i,j} - 1/T_sym * diff(T_sym,x_sym(k,l)) * x_sym(k,l+2) * x_sym(i,j+2) ...
                    + 1/(2*T_sym) * diff(T_sym,x_sym(i,j)) * x_sym(k,l+2)^2;
             end
         end
         F{i,j} = @(x) eval(subs(f{i,j}, x_sym, reshape(x,[4,3])'));
    end
end

start = 0;
terminal = 10;
div_N = 50;
count_div = 10;
plot(init(1), init(2));
hold on;
for count = 1:count_div
    fprintf("count num : %d\n" , count);
    start_p = start / count_div * count;
    terminal_p = terminal / count_div    * count;

    t = linspace(start_p, terminal_p, terminal_p*div_N);
    [t, x] = RK4(@g, t, init); 
    init = x(length(x),:);
    plot(x(:,1),x(:,2),'linewidth',2);
    plot(x(:,5),x(:,6),'linewidth',2);
    plot(x(:,9),x(:,10),'linewidth',2);
end
hold off;

function dxdt = g(t, x) 
    global F; 
    dxdt = zeros(size(x));
    for i = 1:3
        for j = 1:2
             dxdt((i-1)*4+j) = x((i-1)*4+j+2);
             dxdt((i-1)*4+j+2) = F{i,j}(x);
        end
    end
end

function [t, X] = RK4(func, t, prior_x)
    dt = t(2) - t(1);
    nt = length(t);
    X = zeros([nt, numel(prior_x)]);
    X(1,:) = prior_x;
    for i = 1:nt-1
        k1 = func(t(i), X(i,:)); 
        k2 = func(t(i) + dt/2, X(i,:) + k1*dt/2); 
        k3 = func(t(i) + dt/2, X(i,:) + k2*dt/2); 
        k4 = func(t(i) + dt, X(i,:) + k3*dt); 
        X(i+1,:) = X(i,:) + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
    end
end

