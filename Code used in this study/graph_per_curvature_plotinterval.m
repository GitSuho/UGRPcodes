%%% 1-body condition %%%
close all
clc
hold on;



% numerical discrete points
% load_file = "PPCP_RK1_09261106.txt";
% load_file = "PPCP_RK4_10011044.txt";
% load_file = "PPCP_RK4_11020011.txt";
% load_file = "PPCP_RK2_11020044.txt";

% load_file = "PPEP_RK2_11021318.txt";
% load_file = "PPEP_RK2_11021354.txt";
load_file  = "PPEP_RK4_11091106.txt";

% curv_lis = [0.001:0.05:2.001];
% plit_lis = [0.001:0.05:2.001];

curv_lis = [0.001:0.05:4.501];
plit_lis = [0.001:0.005:0.2501];

nume_data = load(load_file); % column - plot interval, row - curvature 

% %for PPCP
% for i = 1:5:length(curv_lis)
%     for j = 1:5:length(plit_lis)
%         scatter3(curv_lis(i), plit_lis(j), nume_data(i, j), 'filled', 'ro');  
%     end
% end

%for PPEP
for i = 1:length(nume_data)
    scatter3(nume_data(i, 1), nume_data(i, 2), nume_data(i, 3), 'filled', 'ro'); 
end


xlabel("curvature (rad/m)");
ylabel("plot interval (m)");
zlabel("relative error (%)");
title("Relative Error of curvature and plot interval")



% %theoritical graph
% curv_info = [0, 2, 300];
% plti_info = [0, 2, 300];
% 
% [curv , plti] = meshgrid(curv_info(1):(curv_info(2)-curv_info(1))/curv_info(3):curv_info(2) , ...
%                          plti_info(1):(plti_info(2)-plti_info(1))/plti_info(3):plti_info(2));
% re_er = NumeErr_EU1(curv, plti);
% mesh(curv, plti, re_er);


for i = 1:5:length(curv_lis)
    for j = 1:5:length(plit_lis)
        scatter3(curv_lis(i), plit_lis(j), NumeErr(curv_lis(i), plit_lis(j)), 'filled', 'bo');
    end
end



hold off;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function rela_err = NumeErr_EU1(curvature, plot_interval)
    %2nd order
    % r = 1./curvature; d = plot_interval;
    % simp_err = sqrt((r - d.^2./sqrt(4.*r.^2 + d.^2)).^2 + ( 2.*r.*d./sqrt(4.*r.^2 + d.^2) ).^2 ) - r;
    
    %1st order
    simp_err = sqrt(plot_interval.^2 + 1./curvature.^2) - 1./curvature;
    rela_err = 100.*simp_err.*curvature;
end

%Theoritical error function of the numerical method
function rela_err = NumeErr(curvature, plot_interval)
    % simp_err = sqrt(plot_interval.^2 + 1./curvature.^2) - 1./curvature;
    % rela_err = 100.*simp_err.*curvature;
    
    r = 1/ curvature;

    d = plot_interval;

    k1 = [0, d];

    x2 = [r - k1(1)/2 , k1(2)/2];
    dx2 = d*x2(2)/sqrt(x2(1).^2 + x2(2).^2);
    dy2 = d*x2(1)/sqrt(x2(1).^2 + x2(2).^2);
    k2 = [dx2, dy2];

    x3 = [r - k2(1)/2 , k2(2)/2];
    dx3 = d*x3(2)/sqrt(x3(1).^2 + x3(2).^2);
    dy3 = d*x3(1)/sqrt(x3(1).^2 + x3(2).^2);
    k3 = [dx3, dy3];

    x4 = [r - k3(1) , k3(2)];
    dx4 = d*x4(2)/sqrt(x4(1).^2 + x4(2).^2);
    dy4 = d*x4(1)/sqrt(x4(1).^2 + x4(2).^2);
    k4 = [dx4, dy4];

    simp_err = abs( sqrt( ((r - (k1(1) + 2* k2(1) + 2*k3(1) + k4(1))/6)  ).^2  ...
                +  (((k1(2) + 2* k2(2) + 2*k3(2) + k4(2))/6) ).^2) - r );

    rela_err = 100.*simp_err.*curvature;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%