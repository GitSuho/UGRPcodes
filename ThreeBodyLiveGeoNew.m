global m; global G;
clear all
close all
clc

% %Ana Hudomal 16p
% %butterfly 1
% % dx0 = 0.30689; dy0 = 1.2551;
% %butterfly 2
% % dx0 = 0.39295; dy0 = 0.09758;
% %bumblebee
% % dx0 = 0.18428; dy0 = 0.58719;
% x_init = [-1 0;1 0; 0 0]; 
% v_init = [dx0 dy0; dx0 dy0; -2*dx0 -2*dy0]; 
% m = [1 1 1];

% %Xiaoming Lia , Shijun Liao. December 2017 
% I.A_{100}^{i.c.}
% v_1 = 0.0670760777;
% v_2 = 0.5889627892;
%I.A_{77}^{i.c.}
% v_1 = 0.4159559963;
% v_2 = 0.2988672319;
%I.B_{102}^{i.c.}
% v_1 = 0.4784306757;
% v_2 = 0.3771895698;
%II.C_{156}^{i.c.}
% v_1 = 0.3231926176;
% v_2 = 0.3279135713;
% 
% x_init = [-1 0 ; 1 0; 0 0];
% v_init = [v_1 v_2; v_1 v_2; -2*v_1 -2*v_2];
% m = [1, 1, 1];

% %lagrange equliteral triangular solution
% m = [1, 2, 3];
% s3_size = 1; s3_rad = 0;
% x_init = 2*[0 2 ; -sqrt(3) -1 ; sqrt(3) -1];
% v_init = 1/sum(m) * [(-m(2)-m(3)/2) (-sqrt(3)/2*m(3));(m(1)+m(3)/2) (-sqrt(3)/2*m(3)); (m(1)/2-m(2)/2) (sqrt(3)/2*m(1)+sqrt(3)/2*m(2))];
% v_init = s3_size*v_init*[cos(s3_rad) sin(s3_rad) ; -sin(s3_rad) cos(s3_rad)];

%Euler collinear solution
m = [1 2 3];
syms lamb_sym
assume(lamb_sym, 'real')

ecs_equ = (m(1)+m(2))*lamb_sym^5 + (3*m(1)+2*m(2))*lamb_sym^4+(3*m(1)+m(2))*lamb_sym^3-(m(2)+3*m(3))*lamb_sym^2-(2*m(2)+3*m(3))*lamb_sym-(m(2)+m(3)) == 0;
lamb_ecs = solve(ecs_equ, lamb_sym);

fprintf("%f\n", lamb_ecs);

%(m(2)+m(3)*(1+lamb_sym))/(m(1)-m(3)*lamb_sym) - (m(2)-m(3)*(1/(1+lamb_sym))^2)/(m(1)-m(3)*(1/lamb_sym)^2) == 0;
%(-m(1)*m(3)-m(2)*m(3))*lamb_sym^5+(-m(1)*m(3)-2*m(1)*m(3)-2*m(2)*m(3))*lamb_sym^4+(-3*m(1)*m(3)+2*m(3).^2-m(2)*m(3))*lamb_sym^3+(-2*m(1)*m(3)+m(2)*m(3)+3*m(3).^2)*lamb_sym^2+(2*m(2)*m(3)+3*m(3).^2)*lamb_sym+(m(2)*m(3)+m(3).^2) == 0;
%(m(1)+m(2))*lamb_sym^5 + (3*m(1)+2*m(2))*lamb_sym^4+(3*m(1)+m(2))*lamb_sym^3-(m(2)+3*m(3))*lamb_sym^2-(2*m(2)+3*m(3))*lamb_sym-(m(2)+m(3)) == 0;
%(m(2)+m(3)*(1+lamb_sym))/(m(1)-m(3)*lamb_sym) - (m(2)+m(3)*(1/(1+lamb_sym))^2)/(m(1)-m(3)*(1/lamb_sym)^2) == 0;
s_3 = [1 0];
sv_3 = [-1/2 1];
form_ecs = 1/sum(m)*[-m(2)-m(3)-lamb_ecs*m(3) ; m(1)-lamb_ecs*m(3) ; m(1)+lamb_ecs*m(1)+lamb_ecs*m(2)];
x_init = form_ecs*s_3;
v_init = form_ecs*sv_3;


G = 1;
init = [x_init v_init]';
init = init(:)';
hold on;
for i = 1:3
    plot(x_init(i,1), x_init(i,2), 'ko');
end
hold off;

E = 0;
for i = 0:2
    j = mod(i,3)+1;
    k = mod(i+1,3)+1;
    E = E + 1/2*m(j)*norm(v_init(j,:))^2 ...
        - G*m(j)*m(k)/norm(x_init(j,:)-x_init(k,:));
end

syms('x_sym', [3 4]);
assume(x_sym, 'real')

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

global GE;
GE{3, 2} = [];
ge{3, 2} = [];
for i = 1:3
    for j = 1:2
        ge{i, j} = 0;
        for k = 1:3
            for l = 1:2
                ge{i, j} = ge{i, j}  - diff(T_sym, x_sym(k,l))*x_sym(k,l+2)*x_sym(i,j+2) ...
                    + 1/(2*T_sym)*diff(T_sym, x_sym(i,j))*x_sym(k,l+2)^2;
            end
        end
        GE{i,j} = @(x) eval(subs(ge{i,j}, x_sym, reshape(x,[4,3])'));
    end
end



global GEt;
GEt{3,2} = [];
get{3,2} = [];
for i = 1:3
    for j = 1:2
         get{i,j} = T_sym_diff / T_sym * x_sym(i,j+2);
         for k = 1:3
             for l = 1:2
                get{i,j} = get{i,j} - 1/T_sym * diff(T_sym,x_sym(k,l)) * x_sym(k,l+2) * x_sym(i,j+2) ...
                    + 1/(2*T_sym) * diff(T_sym,x_sym(i,j)) * x_sym(k,l+2)^2;
             end
         end
         GEt{i,j} = @(x) eval(subs(get{i,j}, x_sym, reshape(x,[4,3])'));
    end
end

global EL;
EL{3, 2} = [];
el{3, 2} = [];
for i = 1:3 %particles' index
    next_1 = mod(i, 3)+1; 
    next_2 = mod(i+1, 3)+1;
    r1_sq = (x_sym(i, 1) - x_sym(next_1, 1)).^2+(x_sym(i, 2)-(x_sym(next_1, 2))).^2;
    r2_sq = (x_sym(i, 1) - x_sym(next_2, 1)).^2+(x_sym(i, 2)-(x_sym(next_2, 2))).^2;
    for j = 1:2 % x, y 
        el{i, j} = - G*m(next_1)*(x_sym(i, j) - x_sym(next_1, j))/(r1_sq.^(1.5))...
                   - G*m(next_2)*(x_sym(i, j) - x_sym(next_2, j))/(r2_sq.^(1.5));
        EL{i, j} = @(x) eval(subs(el{i,j}, x_sym, reshape(x, [4, 3])'));
    end
end


dt = 1/20;
cycle = 10;
% clear("X")
clear("Xt")
% clear("P")

% X(1, :) = init;
Xt(1, :) = init;
% P(1, :) = init;


hold on;
while(1)
    for i = 1:cycle
        % X(i+1, :) = RK4(X(i, :), @g, dt);
        Xt(i+1, :) = RK4(Xt(i, :), @g_t, dt);
        % P(i+1, :) = RK4(P(i, :), @eulagrange, dt);
    end

    % plot(X(:, 1),X(:, 2), 'r-');
    % plot(X(:, 5),X(:, 6), 'g-');
    % plot(X(:, 9),X(:, 10), 'b-');

    plot(Xt(:, 1),Xt(:, 2), 'y-');
    plot(Xt(:, 5),Xt(:, 6), 'm-');
    plot(Xt(:, 9),Xt(:, 10), 'c-');

    % plot(P(:, 1),P(:, 2), 'k:');
    % plot(P(:, 5),P(:, 6), 'k:');
    % plot(P(:, 9),P(:, 10), 'k:');
    
    % plot_curr_triangle(X(cycle+1,:));
    plot_curr_triangle(Xt(cycle+1,:));
    % plot_curr_triangle(P(cycle+1,:));

    % X(1, :) = X(cycle+1, :);
    Xt(1, :) = Xt(cycle+1, :);
    % P(1, :) = P(cycle+1, :);
    pause(0.1);

end

function null = plot_curr_triangle(pos)
    plot([pos(1), pos(5), pos(9), pos(1)] , [pos(2),pos(6),pos(10),pos(2)] , 'b-');
end

function dxdt = g(x) 
    global GE; 
    dxdt = zeros(size(x));
    for i = 1:3
        for j = 1:2
             dxdt((i-1)*4+j) = x((i-1)*4+j+2);
             dxdt((i-1)*4+j+2) = GE{i,j}(x);
        end
    end
end

function dxdt_t = g_t(x) 
    global GEt; 
    dxdt_t = zeros(size(x));
    for i = 1:3
        for j = 1:2
             dxdt_t((i-1)*4+j) = x((i-1)*4+j+2);
             dxdt_t((i-1)*4+j+2) = GEt{i,j}(x);
        end
    end
end

function slope = eulagrange(x)
    global EL;
    slope = zeros(size(x));
    for i = 1:3
        for j = 1:2
            slope((i-1)*4+j) = x((i-1)*4+j+2);
            slope((i-1)*4+j+2) = EL{i,j}(x);
        end
    end
end

function next_pos = RK4(curr_pos, func, dt)
    k1 = func(curr_pos(:));
    k2 = func(curr_pos(:) + k1*dt/2); 
    k3 = func(curr_pos(:) + k2*dt/2); 
    k4 = func(curr_pos(:) + k3*dt); 
    next_pos = curr_pos(:) + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
end



