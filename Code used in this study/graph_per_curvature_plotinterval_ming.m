%%% 1-body condition %%%
close all
clc
hold on;



% numerical discrete points
load_file = "PPCP_RK4_10011044.txt";
curv_lis = [0.001:0.05:2.001];
plit_lis = [0.001:0.05:2.001];
nume_data = load(load_file); % column - plot interval, row - curvature 

for i = 1:5:length(curv_lis)/2
    for j = 1:5:length(plit_lis)/2
        scatter3(curv_lis(i), plit_lis(j), nume_data(i, j), 'filled', 'ro');
    end
end

xlabel("curvature (rad/m)");
ylabel("plot interval (m)");
zlabel("relative error (%)");
title("Relative Error of curvature and plot interval")


%theoritical graph
curv_info = [0, 1, 100];
plti_info = [0, 1, 100];

[curv , plti] = meshgrid(curv_info(1):(curv_info(2)-curv_info(1))/curv_info(3):curv_info(2) , ...
                         plti_info(1):(plti_info(2)-plti_info(1))/plti_info(3):plti_info(2));


re_er = zeros(101);

for i = 1: 101
    for j = 1: 101
        re_er(i,j) = real(NumeErr(curv(1,i),plti(j,1)));
    end
end

% re_er = NumeErr(curv, plti);
surf(curv, plti, re_er);

hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Theoritical error function of the numerical method
function rela_err = NumeErr(curvature, plot_interval)
    % simp_err = sqrt(plot_interval.^2 + 1./curvature.^2) - 1./curvature;
    % rela_err = 100.*simp_err.*curvature;
    
    r = 1./curvature;
    d = plot_interval;

    m2 = -2.*r./d;

    dx3 = -d./(2.*sqrt(1+m2.^2));
    m3 = -(r + dx3)./(m2.*dx3);

    dx4 = -d./sqrt(1+m3.^2);
    m4 = -(r + dx4)./(m3.*dx4);

    tx2 = -d./sqrt(1+m2.^2);
    tx3 = -d./sqrt(1+m3.^2);
    tx4 = -d./sqrt(1+m4.^2);

    x2 = r + tx2./3 + tx3./3 + tx4./6;
    y2 = d./6 + m2*tx2./3 + m3*tx3./3 + m4.*tx4./6;

    simp_err = sqrt(x2.^2 + y2.^2) - r;
    rela_err = 100.*simp_err.*curvature;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%