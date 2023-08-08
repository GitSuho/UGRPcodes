
%Lagrange's Equilatetal Triangular Solutions
% x_init = [-0.0833 0.7217 ; -0.5833 -0.1443 ; 0.4167 -0.1443];
% v_init = [-2.7678 0.7959 ; -0.6464 -0.4289 ; -0.6464 2.0206];

x_init = [-1 -1;0 0; 1 1]; 
v_init = [0 1; 0 0; 0 -1]; 
global m;
m = [1 1 1];
global G;
G = 10;
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

global EL;
EL{3, 2} = [];
el{3, 2} = [];
for i = 1:3 %particles' index

    for j = 1:2 % x, y 
        next_1 = mod(i, 3)+1; next_2 = mod(i+1, 3)+1;
        r = (x_sym(next_1, 1) - x_sym(i, 1)).^2+(x_sym(next_1, 2)-(x_sym(i, 2)).^2)...
            + (x_sym(next_2, 1) - x_sym(i, 1)).^2+(x_sym(next_2, 2)-(x_sym(i, 2)).^2);
        el{i, j} = G*m(i)*m(next_1)*(x_sym(next_1, j) - x_sym(i, j))/(r.^(1.5))...
                 + G*m(i)*m(next_2)*(x_sym(next_2, j) - x_sym(i, j))/(r.^(1.5));
        EL{i, j} = @(x) eval(subs(el{i,j}, x_sym, reshape(x, [4, 3])'));
    end
end


dt = 1/100;
i = 1;
clear("X")
clear("P")
cycle = 20;
% X = zeros(cycle, 12);
% P = zeros(cycle, 12);

X(1, :) = init;
P(1, :) = init;

% X(cycle, :) = init;
% P(cycle, :) = init;


while(1)
    % X(1,:) = X(cycle,:);
    % P(1,:) = P(cycle,:);
    for foo = 2:cycle
        k1 = g(X(i, :));
        k2 = g(X(i, :) + k1*dt/2); 
        k3 = g(X(i, :) + k2*dt/2); 
        k4 = g(X(i, :) + k3*dt); 
        X(i+1, :) = X(i, :) + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
        
        % l1 = eulagrange(P(i, :));
        % l2 = eulagrange(P(i, :) + k1*dt/2); 
        % l3 = eulagrange(P(i, :) + k2*dt/2); 
        % l4 = eulagrange(P(i, :) + k3*dt); 
        % P(i+1, :) = P(i, :) + (l1 + 2*l2 + 2*l3 + l4)*dt/6;
        i = i + 1;
        
    end
    hold off;
    hold on;
    plot(X(:, 1),X(:, 2), 'r-');
    plot(X(:, 5),X(:, 6), 'g-');
    plot(X(:, 9),X(:, 10), 'b-');
    % plot([X(i, 1), X(i, 5), X(i, 9), X(i, 1)],[X(i, 2), X(i, 6), X(i, 10), X(i, 2)], 'black:');
    plot(P(:, 1),P(:, 2), 'c:', 'linewidth', 2);
    plot(P(:, 5),P(:, 6), 'm:', 'linewidth', 2);
    plot(P(:, 9),P(:, 10), 'y:', 'linewidth', 2);
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




