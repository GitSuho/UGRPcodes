%%% 1-body condition %%%
close all
clc

%Initial consitions
E_const = 0.5;
global ecc;
ecc = 0.5;

[E_1, theta] = meshgrid(-0.2:0.005:0,0:pi/50:2*pi);
r = OrbEqu_r(theta, E_1);

hold on;

mesh(r.*cos(theta), r.*sin(theta), E_1);
% surface(r.*cos(theta), r.*sin(theta), E_1);
% plot3(0, 0,0, 'blacko');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%calculate the r which is interval of two particles
function r = OrbEqu_r(theta, E_const)
    global ecc;
    M = 1;m = 1 ;G = 1;
    v_init = sqrt(E_const./(-m/2+m/(ecc+1)));
    x_init = (ecc+1).*G.*M./(v_init.^2);
    k_const = G*M;
    l_const = x_init.*v_init.*sin(acos((x_init.^2+v_init.^2-(x_init.^2+v_init.^2)) ...
        ./(2.*x_init.*v_init)));%r*v*sin(theta_0)
    A_coeff =  (1./x_init - (k_const./(l_const.^2)));
    r = 1./(A_coeff.*cos(theta)+(k_const./(l_const.^2)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%