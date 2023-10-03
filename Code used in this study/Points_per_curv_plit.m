%%% 1-body condition %%%
clear all;

curv_lis = [0.001:0.05:2.001];
plit_lis = [0.001:0.05:2.001];

file_name = "PPCP_RK4_10011044.txt";
wfile = fopen(file_name, 'w');
ord = 4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for curv = curv_lis
    for plot_interval = plit_lis
        M = 0.1;m = 1 ;G = 1;
        v_init = [0, 1];
        x_init = [1/curv, 0];
        
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
        %Orbital equation's constants
        global l_const; global A_coeff; global k_const;
        k_const = G*M;
        l_const = norm(x_init)*norm(v_init)*sin(acos((norm(x_init).^2+norm(v_init).^2-norm(x_init-v_init).^2) ...
            /(2*norm(x_init)*norm(v_init))));%r*v*sin(theta_0)
        A_coeff =  (1/norm(x_init) - (k_const/(l_const.^2)))/(x_init(1)/norm(x_init));
      
                
        dt = 1/50;       
        dt = Fit_dt(dt, plot_interval ,init, @geo, ord);
        X = RKn(init, @geo, dt, ord);
        X_error = Err_est(X(1), X(2), 1/curv);

        fprintf(wfile ,'%f    ', X_error);
        
    end
    fprintf(wfile , '\n');
end
fclose(wfile);
fprintf('program end\n');

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
    r = norm(OrbEqu(theta));
end
%calculate a degree when x-y coordinate given
function degree = Find_degree(x, y)
    degree = mod (atan2(y, x), 2*pi);
end
%find a dt value that make given plot interval
function result_dt = Fit_dt(dt, fit_val, curr_pos, func, ord)
    accuracy = 0.00000001;
    next_pos = RKn(curr_pos,func, dt, ord);
    disp  = sqrt((next_pos(1) - curr_pos(1)).^2 + (next_pos(2) - curr_pos(2)).^2);
    while( abs(disp - fit_val) > accuracy)
        dt = fit_val*dt/disp;
        next_pos = RKn(curr_pos,func, dt, ord);
        disp  = sqrt((next_pos(1) - curr_pos(1)).^2 + (next_pos(2) - curr_pos(2)).^2);
    end
    result_dt = dt;
end
%calculate error when x-y coordinate is given
function result_err = Err_est(x, y, r_0)
    % rr = OrbEqu_r(Find_degree(x, y));
    result_err = 100 * abs(r_0 - sqrt(x.^2 + y.^2)) / r_0; %relative_err
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%