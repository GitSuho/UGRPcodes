%%% 1-body condition %%%
close all
clc
hold on;



% numerical discrete points
load_file = "PPCP_RK1_09251153.txt";
curv_lis = [0.001:0.05:2.001];
plit_lis = [0.001:0.05:2.001];
nume_data = load(load_file); % column - plot interval, row - curvature 

for i = 1:5:length(curv_lis)
    for j = 1:5:length(plit_lis)
        scatter3(plit_lis(j), curv_lis(i), nume_data(i, j), 'filled', 'ro');
    end
end




%theoritical graph
curv_info = [0, 2, 100];
plti_info = [0, 2, 100];

[curv , plti] = meshgrid(curv_info(1):(curv_info(2)-curv_info(1))/curv_info(3):curv_info(2) , ...
                         plti_info(1):(plti_info(2)-plti_info(1))/plti_info(3):plti_info(2));
re_er = NumeErr(curv, plti);
mesh(curv, plti, re_er);

hold off;
xlabel("curvature (rad/m)");
ylabel("plot interval (m)");
zlabel("relative error (%)");
title("Relative Error of curvature and plot interval")




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Theoritical error function of the numerical method
function rela_err = NumeErr(curvature, plot_interval)
    simp_err = sqrt(plot_interval.^2 + 1./curvature.^2) - 1./curvature;
    rela_err = 100.*simp_err.*curvature;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%