global m; global G;
clear all
close all
clc

%lagrange equliteral triangular solution
m = [1, 2, 3];
s3_size = 1; s3_rad = 0;
x_init = 2*[0 2 ; -sqrt(3) -1 ; sqrt(3) -1];
v_init = 1/sum(m) * [(-m(2)-m(3)/2) (-sqrt(3)/2*m(3));(m(1)+m(3)/2) (-sqrt(3)/2*m(3)); (m(1)/2-m(2)/2) (sqrt(3)/2*m(1)+sqrt(3)/2*m(2))];
v_init = s3_size*v_init*[cos(s3_rad) sin(s3_rad) ; -sin(s3_rad) cos(s3_rad)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%basic setting
G = 1;
init = [x_init v_init]';
init = init(:)';
hold on;
for i = 1:3
    plot(x_init(i,1), x_init(i,2), 'ko');
end
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

%Geodesic solution which domain is arclength
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
%Euler-Lagrange solution
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


%%%%
T_sym_diff = 0;
for i = 0:2
    j = mod(i,3)+1;
    k = mod(i+1,3)+1;
    T_sym_diff = T_sym_diff + G * m(j) * m(k) * dot(x_sym(j,1:2)-x_sym(k,1:2),x_sym(j,3:4)-x_sym(k,3:4)) ...
                                / norm(x_sym(j,1:2)-x_sym(k,1:2))^2 ;
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
%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_name = "test_threebody.txt";
wfile = fopen(file_name, 'w');
fprintf(wfile,"%s\n", file_name);
fprintf(wfile, "initial value :\n");
for i = 0:2
    fprintf(wfile, "[%f , %f , %f , %f]\n", init(i*4+1), init(i*4+2), init(i*4+3), init(i*4+4));
end

dt = 1/10;
cycle = 75;

clear("X")
clear("P")
X(1, :) = init;
P(1, :) = init;

%%%%
clear("Xt")
Xt(1, :) = init;
%%%

hold on;
count = 0;
for var = 1:20
    fprintf("count %d\n", var);
    %expect next position used by numerical method and calculate error value
    for i = 1:cycle
        count = count + 1;
        fprintf(wfile, "%4d | ", count);
        X(i+1, :) = RK4(X(i, :), @g, dt);
        fprintf(wfile, "Geod_err : %f | ", Err_est_tri(X(i,:)));
        P(i+1, :) = RK4(P(i, :), @eulagrange, dt);
        fprintf(wfile, "Newt_err : %f \n", Err_est_tri(P(i,:)));

        %%%%
        Xt(i+1, :) = RK4(Xt(i, :), @g_t, dt);
        fprintf(wfile, "t_Ge_err : %f \n", Err_est_tri(Xt(i,:)));
        %%%%

    end
    
    %plot points and draw triangles according to last points
    plot(X(:, 1),X(:, 2), 'r--o');
    plot(X(:, 5),X(:, 6), 'g--o');
    plot(X(:, 9),X(:, 10), 'b--o');
    plot(P(:, 1),P(:, 2), 'm--^');
    plot(P(:, 5),P(:, 6), 'c--^');
    plot(P(:, 9),P(:, 10), 'k--^');
    plot_curr_triangle(X(cycle+1,:));
    plot_curr_triangle(P(cycle+1,:));

    %%%%
    plot(Xt(:, 1),Xt(:, 2), 'r*');
    plot(Xt(:, 5),Xt(:, 6), 'r*');
    plot(Xt(:, 9),Xt(:, 10), 'r*');
    plot_curr_triangle(Xt(cycle+1,:));
    Xt(1, :) = Xt(cycle+1, :);
    %%%%



    X(1, :) = X(cycle+1, :);
    P(1, :) = P(cycle+1, :);
    pause(0.1);
end

hold off;
fclose(wfile);
fprintf("program end\n");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate error which points should make equliteral triangular solution
function err_result = Err_est_tri(pos)
    r1 = norm(pos(1:2)-pos(5:6));
    r2 = norm(pos(5:6)-pos(9:10));
    r3 = norm(pos(9:10)-pos(1:2));
    err_result = sum(abs(r1-r2)+abs(r2-r3)+abs(r3-r1))/3;
end
%draw triangle by given points
function plot_curr_triangle(pos)
    plot([pos(1), pos(5), pos(9), pos(1)] , [pos(2),pos(6),pos(10),pos(2)] , 'y:');
end
%geodesic equation
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

%%%%
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
%%%%


%newtonian equation
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
%Runge-Kutta 4th Order
function next_pos = RK4(curr_pos, func, dt)
    k1 = func(curr_pos(:));
    k2 = func(curr_pos(:) + k1*dt/2); 
    k3 = func(curr_pos(:) + k2*dt/2); 
    k4 = func(curr_pos(:) + k3*dt); 
    next_pos = curr_pos(:) + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
end