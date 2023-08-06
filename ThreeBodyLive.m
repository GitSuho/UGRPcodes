
%Lagrange's Equilatetal Triangular Solutions
x_init = [-0.0833 0.7217 ; -0.5833 -0.1443 ; 0.4167 -0.1443];
v_init = [-2.7678 0.7959 ; -0.6464 -0.4289 ; -0.6464 2.0206];

%x_init = [0 2; sqrt(3) -1; -sqrt(3) -1]; 
%v_init = [-sqrt(3.25) 0; sqrt(3) 1; sqrt(3) -1]; 

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

hold on;
dt = 1/div_N;
X(1, :) = init;
i = 1;
while(1)

    k1 = g(X(i, :));
    k2 = g(X(i, :) + k1*dt/2); 
    k3 = g(X(i, :) + k2*dt/2); 
    k4 = g(X(i, :) + k3*dt); 
    X(i+1, :) = X(i, :) + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
    i = i + 1;

    plot(X(:,1), X(:,2));
    plot(X(:,5), X(:,6));
    plot(X(:,9), X(:,10));
    pause(10);

end
hold off;

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
