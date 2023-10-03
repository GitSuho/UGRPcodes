%%% 1-body condition %%%
clear all;


E_list = [-0.05, -0.1, -0.15, -0.20, -0.25];
ecc_list = [0, 0.4, 0.8];
numerical_method = ["EU1","HE2","KU3","RK4"];

for ord = 1:1
    for ii = 2:2
        for jj = 2:2

E = E_list(ii);
ecc = ecc_list(jj);
nume_name = numerical_method(ord);

filename = sprintf("OneBody5_E%1.2f_ecc%1.1f_%s.txt", E, ecc, nume_name);
fprintf(filename);


M = 9;m = 1 ;G = 1;
v_init = [0, sqrt(-E/(-m/2+m/(ecc+1)))];
x_init = [(ecc+1)*G*M/(v_init(2).^2), 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms('x_sym', [2, 2]);
assume(x_sym, 'real');
E = 1/2*m*norm(v_init)^2 - G*M*m/norm(x_init);
T = E + G*M*m/norm(x_sym(1,:));
init = [x_init v_init];
init = init(:);
%Geodesic solution which domain is arclength
global GE;
GE{2} = [];
ge{2} = [];
for i = 1:2
    ge{i} = -(diff(T, x_sym(1,1))*x_sym(2,1) + diff(T, x_sym(1,2))*x_sym(2,2))*x_sym(2,i) ...
       + 1/(2*T)*diff(T, x_sym(1,i))*norm(x_sym(2,:))^2;
    GE{i} = @(x) eval(subs(ge{i}, x_sym, reshape(x,[2,2])'));
end
%Euler-Lagrange solution
global EL;
EL{2} = [];
el{2} = [];
for i = 1:2
    el{i} = -G*M*x_sym(1, i)/(norm(x_sym(1,:))^3);
    EL{i} = @(x) eval(subs(el{i}, x_sym, reshape(x,[2,2])'));
end

global a; global b; global r_0;
r_0 = x_init(1);
a = x_init(1)*(1 + (1+ecc)/(1-ecc))/2;
b = a*sqrt(1-ecc.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt_ind = [0.2, 0.1, 0.05, 0.025];
wfile = fopen(filename, 'w');
fprintf(wfile, "ecc%f,a%f,b%f,r_0,%f\n", ecc, a, b, x_init(1));

for dt_ind = 1:4
    dt = dt_ind(dt_ind);


    X(1, :) = init;
    X_error(1, :) = [0, 0, 0];
    
    
    plot_count = 0;
    while(1)
        for i = 2:101
            X(i, :) = RKn(X(i-1,:), @geo, dt, ord);    
            X_error(i, :) = Err_est(X(i,1), X(i,2));    
    
            if (X_error(i, 1) < X_error(i-1, 1))
                break;
            end
            plot_count = plot_count + 1;
        end
        
        if(X_error(101, 1) < X_error(100, 1) || i < 101)
            break;
        end
    
        X(1, :) = X(101, :);
        X_error(1, :) = X_error(101, :);  
    end

    newt_dt = Newt_stepsize_pred(plot_count);
    fprintf(wfile, "plot_num_perOnePeriod%d,geod_stepsize%f,newt_step_size%f\n", plot_count, dt, newt_dt);

end

fclose(wfile);

        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%geodesic equation
function dxdt = geo(x)
    global GE;
    dxdt = zeros(size(x));
    for i = 1:2
        dxdt(i) = x(i+2);
        dxdt(i+2) = GE{i}(x);
    end
end
%newtonian equation
function slope = lag(x)
    global EL;
    slope = zeros(size(x));
    for i = 1:2
        slope(i) = x(i+2);
        slope(i+2) = EL{i}(x);
    end
end
%Runge-Kutta nth Order
function next_pos = RKn(curr_pos, func, dt, order)
    switch order
        case 4
            k1 = func(curr_pos);
            k2 = func(curr_pos + k1*dt/2); 
            k3 = func(curr_pos + k2*dt/2); 
            k4 = func(curr_pos + k3*dt); 
            next_pos = curr_pos + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
        case 3
            k1 = func(curr_pos);
            k2 = func(curr_pos + k1*dt/2);
            k3 = func(curr_pos - k1*dt + 2*k2*dt);
            next_pos = curr_pos + (k1 + 4*k2 + k3)*dt/6;
        case 2
            k1 = func(curr_pos);
            k2 = func(curr_pos + k1*dt);
            next_pos = curr_pos + (k1+k2)*dt/2;
        case 1
            next_pos = curr_pos + func(curr_pos)*dt;
    end
end
%calculate a degree when x-y coordinate given
function degree = Find_degree(x, y)
    if (y >= 0 )
        degree = acos( x/sqrt(x.^2 + y.^2));
    else
        degree = 2*pi - acos( x/sqrt(x.^2 + y.^2));
    end
end
%calculate error when x-y coordinate is given
function result_err = Err_est(x, y)
    theta = Find_degree(x, y);
    simple_err = 0;%abs(OrbEqu_r(theta) - norm(x, y));
    relative_err = 0;%simple_err / OrbEqu_r(theta) * 100;
    result_err = [theta, simple_err, relative_err];
end
function step_size = Newt_stepsize_pred(plot_num)
    global a ; global b; global G; global M;global r_0;
    step_size = ((2*b*pi)/plot_num)*sqrt((2*a.^3)/(r_0*G*M*(2*a-r_0)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%