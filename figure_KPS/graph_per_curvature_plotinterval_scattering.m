%%% 1-body condition %%%
close all
clc
hold on;



% numerical discrete points
load_file = "PPCP_RK1_102209.txt";
nume_data = load(load_file); % column - plot interval, row - curvature 
sz = 20;


% for j = 1:length(nume_data)/4
%     i = j;
%     scatter3(nume_data(i,2), nume_data(i,1), nume_data(i,3), sz, 'filled', 'MarkerFaceColor' , [1 0 0]);
% end
% for j = 1:length(nume_data)/4
%     i = j + length(nume_data)/4;
%     scatter3(nume_data(i,2), nume_data(i,1), nume_data(i,3), sz, 'filled', 'MarkerFaceColor' , [1 0.5 0]);
% end
% for j = 1:length(nume_data)/4
%     i = j + 2*length(nume_data)/4;
%     scatter3(nume_data(i,2), nume_data(i,1), nume_data(i,3), sz, 'filled', 'MarkerFaceColor' , [1 1 0]);
% end
% for j = 1:length(nume_data)/4
%     i = j + 3*length(nume_data)/4;
%     scatter3(nume_data(i,2), nume_data(i,1), nume_data(i,3), sz, 'filled', 'MarkerFaceColor' , [0 1 0]);
% end


scatter3(nume_data(:,2), nume_data(:,1), nume_data(:,3), sz, 'filled', 'MarkerFaceColor' , [1 0 0],'Displayname', 'e = 0.9');

load_file = "PPCP_RK1_102206.txt";
nume_data = load(load_file); % column - plot interval, row - curvature 


scatter3(nume_data(:,2), nume_data(:,1), nume_data(:,3), sz, 'filled', 'MarkerFaceColor' , [1 0.5 0],'Displayname', 'e = 0.6');

load_file = "PPCP_RK1_102203.txt";
nume_data = load(load_file); % column - plot interval, row - curvature 

scatter3(nume_data(:,2), nume_data(:,1), nume_data(:,3), sz, 'filled', 'MarkerFaceColor' , [1 1 0],'Displayname', 'e = 0.3');

legend;
xlabel("curvature (rad/m)");
ylabel("plot interval (m)");
zlabel("absolute error (m)");
title("Absolute Error of curvature and plot interval")


%theoritical graph
curv_info = [0, 2.5, 50];
plti_info = [0, 2.5, 50];

[curv , plti] = meshgrid(curv_info(1):(curv_info(2)-curv_info(1))/curv_info(3):curv_info(2) , ...
                         plti_info(1):(plti_info(2)-plti_info(1))/plti_info(3):plti_info(2));
re_er = NumeErr(curv, plti);
s = surf(curv, plti, re_er,'FaceAlpha', 0.5,'Displayname','');
s.EdgeColor = 'none';
view (-45, 30);
hold off;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Theoritical error function of the numerical method
function rela_err = NumeErr(curvature, plot_interval)
    simp_err = sqrt(plot_interval.^2 + 1./curvature.^2) - 1./curvature;
%     rela_err = 100.*simp_err.*curvature;
rela_err = simp_err;
    
    % r = 1./curvature;
    % d = plot_interval;
    % 
    % m2 = -2.*r./d;
    % 
    % dx3 = -d./(2.*sqrt(1+m2.^2));
    % m3 = -(r + dx3)./(m2.*dx3);
    % 
    % dx4 = -d./sqrt(1+m3.^2);
    % m4 = -(r + dx4)./(m3.*dx4);
    % 
    % tx2 = -d./sqrt(1+m2.^2);
    % tx3 = -d./sqrt(1+m3.^2);
    % tx4 = -d./sqrt(1+m4.^2);
    % 
    % x2 = r + tx2./3 + tx3./3 + tx4./6;
    % y2 = d./6 + m2*tx2./3 + m3*tx3./3 + m4.*tx4./6;
    % 
    % simp_err = sqrt(x2.^2 + y2.^2) - r;
    % rela_err = 100.*simp_err.*curvature;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%