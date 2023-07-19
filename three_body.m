%% 3-body problem in a geometric view

x_init = [0 0;1 0;1 2]; % [x1 x2; x3 x4; x5 x6] init position
v_init = [1 0;0 1;-1/sqrt(2) -1/sqrt(2)]; % [x1' x2'; x3' x4'; x5' x6'] init velocity
m = [1 1 1]; % masses
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

% [x1  x2  x1' x2']
% [x3  x4  x3' x4']
% [x5  x6  x5' x6']

syms('x_sym', [3 4]);
assume(x_sym, 'real')

% T = E - V
T_sym = E;
for i = 0:2
    j = mod(i,3)+1;
    k = mod(i+1,3)+1;
    T_sym = T_sym + G*m(j)*m(k)/norm(x_sym(j,1:2)-x_sym(k,1:2));
end

F{3,2} = [];
f{3,2} = [];
for i = 1:3
    for j = 1:2
         f{i,j} = 0;
         for k = 1:3
             for l = 1:2
                f{i,j} = f{i,j} - (diff(T_sym,x_sym(k,l)) * x_sym(i,j+2) ...
                    + 1/2 * diff(T_sym,x_sym(i,j)) * x_sym(k,l+2)) * x_sym(k,l+2);
             end
         end
         f{i,j} = f{i,j} / T_sym;
         F{i,j} = @(x) eval(subs(f{i,j}, x_sym, reshape(x,[4,3])'));
    end
end

g = @(t,x) [x(3);
            x(4);
            F{1,1}(x);
            F{1,2}(x);
            x(7);
            x(8);
            F{2,1}(x);
            F{2,2}(x);
            x(11);
            x(12);
            F{3,1}(x);
            F{3,2}(x)];
N = 1000;
t = linspace(0, 30, N);
[~, x] = ode45(g, t, init);

% plot
plot(x(:,1),x(:,2),'linewidth',2);
hold on;
plot(x(:,5),x(:,6),'linewidth',2);
plot(x(:,9),x(:,10),'linewidth',2);
hold off;

%%
for i = 1:N
    plot(x(1:i,1),x(1:i,2),'linewidth',2);
    hold on;
    plot(x(1:i,5),x(1:i,6),'linewidth',2);
    plot(x(1:i,9),x(1:i,10),'linewidth',2);
    hold off;
    drawnow;
end
%%
% x_k'' = 1/T*(T)_(x_i)*(x_i')*(x_k') + 1/(2T)*(T)_(x_k)(x_i')^2

clear all
clc

x_init = [0 0;.1 0;.1 .1]; % [x1 x2; x3 x4; x5 x6] init position
v_init = [3 0;-3 3;0 3]; % [x1' x2'; x3' x4'; x5' x6'] init velocity
m = [1 1 1]; % masses
G = 9.8;
init = [x_init v_init]';
init = init(:)';
init = [0     0     3     0     1     0    -3     3     1     1     0     3];
E = 0;
for i = 0:2
    j = mod(i,3)+1;
    k = mod(i+1,3)+1;
    E = E + 1/2*m(j)*norm(v_init(j,:))^2 ...
        - G*m(j)*m(k)/norm(x_init(j,:)-x_init(k,:));
end
syms('x_sym', [3 4]);
% [x1  x2  x1' x2']
% [x3  x4  x3' x4']
% [x5  x6  x5' x6']
T_sym = E;
for i = 0:2
    j = mod(i,3)+1;
    k = mod(i+1,3)+1;
    T_sym = T_sym + G*m(j)*m(k)/norm(x_sym(j,1:2)-x_sym(k,1:2));
end
% T = E - V

F{3,2} = [];
f{3,2} = [];
for i = 1:3
    for j = 1:2
         f{i,j} = 0;
         for k = 1:3
             for l = 1:2
                f{i,j} = f{i,j} - 1/T_sym * diff(T_sym,x_sym(k,l)) * x_sym(k,l+2) * x_sym(i,j+2) ...
                    + 1/(2*T_sym) * diff(T_sym,x_sym(i,j)) * x_sym(k,l+2)^2;
             end
         end
         F{i,j} = @(x) eval(subs(f{i,j}, x_sym, reshape(x,[4,3])'));
    end
end
diff(T_sym,x_sym(1,2))

for i = 1:3
    for j = 1:2
        eval(subs(diff(T_sym,x_sym(1,2)), x_sym, [x_init v_init]))
    end
end

g = @(t,x) [x(3);
            x(4);
            F{1,1}(x);
            F{1,2}(x);
            x(5);
            x(6);
            F{2,1}(x);
            F{2,2}(x);
            x(9);
            x(10);
            F{3,1}(x);
            F{3,2}(x)];
N = 10000;
t = linspace(0, 2, N);
[t, x] = ode45(g, t, init);
% plot

plot(x(:,1),x(:,2),'linewidth',2);
hold on;
plot(x(:,5),x(:,6),'linewidth',2);
plot(x(:,9),x(:,10),'linewidth',2);
hold off;

