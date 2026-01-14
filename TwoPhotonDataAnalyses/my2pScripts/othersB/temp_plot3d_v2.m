close all


if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
    x = beta_memoryPrecision_summary;
    y = beta_choiceMemory_summary;
    z = beta_choiceMemory_baseline_summary;
    tempStr = 'trial-level';
elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
    x = beta_seqPrecision_summary;
    y = beta_gProb_summary;
    z = beta_gProb_baseline_summary;
    tempStr = 'seq-level';
elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
    x = beta_seqPrecision_summary;
    y = beta_gProb_summary;
    z = beta_choiceMemory_baseline_summary;
    tempStr = 'mix-level';
end


[x_min,x_max] = bounds(x);
[y_min,y_max] = bounds(y);
[z_min,z_max] = bounds(z);


% x_min = -2;
% x_max = 2;
% y_min = 0;
% y_max = 100;
% z_min = -2;
% z_max = 2;

% y=0:0.1:80;
% x = y/50.*cos(y);
% z = y/50.*sin(y);

sampleCount = length(y);

tempSize = 3;%10
scatter3(x,y,z,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
grid on
% xlabel('x');
% ylabel('y');
% zlabel('z');
xlim([x_min x_max]);
ylim([y_min y_max]);
zlim([z_min z_max]);
hold on
% scatter3(x,y_max*ones(1,sampleCount),z,tempSize,'filled','MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7); % project in x-z axis
% hold on
% scatter3(x_max*ones(1,sampleCount),y,z,tempSize,'filled','MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7); % project in y-z axis
% hold on
% scatter3(x,y,z_min*ones(1,sampleCount),tempSize,'filled','MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7); % project in y-z axis
% hold on

set(gca,'linewidth',1.5)
set(gca, 'FontSize', 10)
set(gca,'box','off');% 取消右、上边框
xlabel('Memory precision', 'FontSize', 10, 'FontWeight', 'normal');
ylabel('Meta-memory', 'FontSize', 10, 'FontWeight', 'normal');
zlabel('Baseline', 'FontSize', 10, 'FontWeight', 'normal');



% x_min = -2;
% x_max = 2;
% y_min = 0;
% y_max = 100;
% z_min = -2;
% z_max = 2;
% 
% y=0:0.1:80;
% x = y/50.*cos(y);
% z = y/50.*sin(y);
% 
% sampleCount = size(y,2);
% 
% plot3(x,y,z, 'LineWidth', 2)
% grid on
% xlabel('x');
% ylabel('y');
% zlabel('z');
% xlim([x_min x_max]);
% ylim([y_min y_max]);
% zlim([z_min z_max]);
% hold on
% plot3(x, y_max*ones(1,sampleCount), z, 'LineWidth', 2); % project in x-z axis at y=100
% plot3(x_max*ones(1,sampleCount), y, z, 'LineWidth', 2); % project in y-z axis at x=2
% plot3(x, y, z_min*ones(1,sampleCount), 'LineWidth', 2); % project in y-z axis at z=-2