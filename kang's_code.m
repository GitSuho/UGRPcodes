syms('x', [2, 2]);
%[x1_1 x1_2  = [x1 x2
% x2_1 x2_2]    v1 v2]
assume(x, 'real');
x0 = [1, 0];
v0 = [0, 1];
M = 1; m = 1; G = 9.8;
E = 1/2*m*norm(v0)^2 - G*M*m/norm(x0);
T = E + G*M*m/norm(x(1,:));
F1 = -(diff(T, x(1,1))*x(2,1) + diff(T, x(1,2))*x(2,2))*x(2,1) ...
       + 1/(2*T)*diff(T, x(1,1))*norm(x(2,:))^2;
F2 = -(diff(T, x(1,1))*x(2,1) + diff(T, x(1,2))*x(2,2))*x(2,2) ...
       + 1/(2*T)*diff(T, x(1,2))*norm(x(2,:))^2;
f1 = @(t) eval(subs(F1, x, reshape(t,[2,2])'));
f2 = @(t) eval(subs(F2, x, reshape(t,[2,2])'));
% [x1 x2 v1 v2] -> [x1 x2
%                   v1 v2]
g = @(t,x) [x(3);
            x(4);
            f1(x);
            f2(x)];
N = 1000;
for i = 1 : 1
    [~, x] = ode45(g, linspace(0, 1, N), [x0 v0]);
    plot(x(:,1),x(:,2),'linewidth',3);
    axis equal
    hold on ;
    x0 = x(end,1:2);
    v0 = x(end,3:4)/norm(x(end,3:4));
end
%%
% classical approach
% x = [x1 x2 x1' x2']
x0 = [1, 0];
v0 = [0, 1]*(E+G*M*m/norm(x0));
g = @(t, x) [
    x(3);
    x(4);
    -G*M/norm(x(1:2))^3*x(1);
    -G*M/norm(x(1:2))^3*x(2)];
[t, x] = ode45(g, linspace(0, .2, N), [x0 v0]);
plot(x(:,1),x(:,2),'linewidth',3);
axis([0,1,0,1])
hold off