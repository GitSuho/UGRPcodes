global m; global G;

%Lagrange's Equilatetal Triangular Solutions
% x_init = [-0.0833 0.7217 ; -0.5833 -0.1443 ; 0.4167 -0.1443];
% v_init = [-2.7678 0.7959 ; -0.6464 -0.4289 ; -0.6464 2.0206];

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
% I.Z_{100}^{i.c.}
% v_1 = 0.0670760777;
% v_2 = 0.5889627892;
%I.A_{77}^{i.c.}
% v_1 = 0.4159559963;
% v_2 = 0.2988672319;
%I.B_{102}^{i.c.}
% v_1 = 0.4784306757;
% v_2 = 0.3771895698;
%II.C_{156}^{i.c.}
v_1 = 0.3231926176;
v_2 = 0.3279135713;

x_init = [-1 0 ; 1 0; 0 0];
v_init = [v_1 v_2; v_1 v_2; -2*v_1 -2*v_2];
m = [1, 1, 1];


% x_init = [2 1 ; 2 2; 3 2];
% v_init = [-1 0 ; -sqrt(2) sqrt(2);1 0 ];
% m = [1, 1, 1];


G = 1;
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


dt = 1/100;
 cycle = 50;
clear("X")
clear("P")

X(1, :) = init;
% P(1, :) = init;


while(1)
    for i = 2:cycle
        X(i+1, :) = RK4(X(i, :), @g, dt);
        % P(i+1, :) = RK4(P(i, :), @eulagrange, dt);
    end
    hold off;
    hold on;
    plot(X(:, 1),X(:, 2), 'b-');
    plot(X(:, 5),X(:, 6), 'r-');
    plot(X(:, 9),X(:, 10), 'black-');

    % plot(P(:, 1),P(:, 2), 'c:', 'linewidth', 2);
    % plot(P(:, 5),P(:, 6), 'm:', 'linewidth', 2);
    % plot(P(:, 9),P(:, 10), 'y:', 'linewidth', 2);
    X(1, :) = X(cycle+1, :);
    % P(1, :) = P(cycle+1, :);
    pause(0.1);

end

function dxdt = g(x) 
    global F; 
    dxdt = zeros(size(x));
    for i = 1:3
        for j = 1:2
             dxdt((i-1)*4+j) = x((i-1)*4+j+2);
             dxdt((i-1)*4+j+2) = F{i,j}(x);
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



