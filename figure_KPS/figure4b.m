%%% 1-body condition %%%
clear all;



plot_interval = 0.5;
cont_curvature = 1/5;

E = -2*cont_curvature;
ord = 4;
M = 1;m = 1 ;G = 1;
v_init = [0, sqrt(-E*2)];
x_init = [1/(-E*2), 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms('x_sym', [2, 2]);
assume(x_sym, 'real');
E = 1/2*m*norm(v_init)^2 - G*M*m/norm(x_init);
T = E + G*M*m/norm(x_sym(1,:));
init1 = [x_init v_init];
init1 = init1(:);

init2 = [x_init -v_init];
init2 = init2(:);

%Geodesic solution which domain is arclength
global GE;
GE{2} = [];
ge{2} = [];
for i = 1:2
    ge{i} = -(diff(T, x_sym(1,1))*x_sym(2,1) + diff(T, x_sym(1,2))*x_sym(2,2))*x_sym(2,i) ...
       + 1/(2*T)*diff(T, x_sym(1,i))*norm(x_sym(2,:))^2;
    GE{i} = @(x) eval(subs(ge{i}, x_sym, reshape(x,[2,2])'));
end
%Orbital equation's constants
global l_const; global A_coeff; global k_const;
k_const = G*M;
l_const = norm(x_init)*norm(v_init)*sin(acos((norm(x_init).^2+norm(v_init).^2-norm(x_init-v_init).^2) ...
    /(2*norm(x_init)*norm(v_init))));%r*v*sin(theta_0)
A_coeff =  (1/norm(x_init) - (k_const/(l_const.^2)))/(x_init(1)/norm(x_init));

global a_const; global b_const;
a_const = x_init(1)*(1 + (1+0)/(1-0))/2;
b_const = a_const*sqrt(1-0.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dt = plot_interval / v_init(2);
X2 = RKn(init1, @geo, dt, ord);
% X2 = RKn(X2_0, @geo, dt, ord);
X0 = RKn(init2, @geo, dt, ord);
% X0 = RKn(X0_0, @geo, dt, ord);

disc_curvature = Find_discurv(X2(1:2) ,x_init, X0(1:2), plot_interval);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find discrete curvature
function curv = Find_discurv(X2, X1, X0, d)
X_mid = [(X2(1)+X0(1))/2 (X2(2)+X0(2))/2];    
aa = norm(X_mid - X1);
bb = norm(X1 - X0);
cc = norm(X_mid - X0);
theta_inter= acos((aa.^2 + bb.^2 - cc.^2) / (2*aa*bb));

curv = 2*cos(theta_inter)/d;
end
%geodesic equation
function dxdt = geo(x)
    global GE;
    dxdt = zeros(size(x));
    for i = 1:2
        dxdt(i) = x(i+2);
        dxdt(i+2) = GE{i}(x);
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
%calculate point of a moving particle when theta is given
function xy_coor = OrbEqu(theta)
    global l_const; global A_coeff; global k_const;
    r = 1/(A_coeff*cos(theta)+(k_const/(l_const.^2)));
    xy_coor = [r*cos(theta), r*sin(theta)];
end
%calculate the r which is interval of two particles
function r = OrbEqu_r(theta)
    global l_const; global A_coeff; global k_const;
    r = 1/(A_coeff*cos(theta)+(k_const/(l_const.^2)));
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%