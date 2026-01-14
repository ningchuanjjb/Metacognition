close all;

x = beta_seqPrecision_summary.^2;
y = beta_gProb_summary.^2;
z = beta_choiceMemory_baseline_summary.^2;


% x = abs(beta_seqPrecision_summary);
% y = abs(beta_gProb_summary);
% z = abs(beta_choiceMemory_baseline_summary);

temp_vectors = [x,y,z];


x = r2_seqPrecision_summary;
y = r2_gProb_summary;
z = r2_choiceMemory_baseline_summary;
% temp_threshold = -1.000;%0,0.01
temp_threshold = 0.1;%0,0.01
temp_vectors_r2 = [x,y,z];

% temp_vectors = [x,y,z];

temp_vectors_n11n = temp_vectors;

temptempBoolIndex = sum(temp_vectors_r2<temp_threshold,2)==3;

% temp_vectors_n11n(temptempBoolIndex,:) = 0;
% temp_vectors_n11n(~temptempBoolIndex,:) = 0;

temp_vectors_n11n = temp_vectors_n11n ./ sum(temp_vectors_n11n,2);

temp_vectors_n11n_valid = temp_vectors_n11n(~isnan(temp_vectors_n11n(:,1)),:);

temp_data = temp_vectors_n11n_valid;
% ternplot(temp_data(:,1),temp_data(:,2),temp_data(:,3));

NUM_AXES_STEPS = 10;
NUM_COLOR_CLASSES = 10;
% figure;
ternplot_pro(temp_data(:,1),temp_data(:,2),temp_data(:,3),NUM_AXES_STEPS,NUM_COLOR_CLASSES);
hold on
% ternplot(temp_data(:,1),temp_data(:,2),temp_data(:,3),'.','color',[1 1 1]*0.5);

hold on
% ternlabel('Memory precision', 'Meta-memory', 'Prior');
% A_label = 'Memory precision';
% B_label = 'Meta-memory';
% C_label = 'Prior';
% text(0.5, -0.05-0.02, A_label, 'horizontalalignment', 'center');
% text(1-0.45*sin(deg2rad(30))+0.02, 0.5+0.02, B_label, 'rotation', -60, 'horizontalalignment', 'center');
% text(0.45*sin(deg2rad(30))-0.02, 0.5+0.02, C_label, 'rotation', 60, 'horizontalalignment', 'center');

% figure;
% Z = ones(size(temp_data,1),1);
% % terncontour(temp_data(:,1),temp_data(:,2),temp_data(:,3), Z);
% % ternpcolor(temp_data(:,1),temp_data(:,2),temp_data(:,3), Z);
% % ternpcolor(temp_data(1:3,1),temp_data(1:3,2),temp_data(1:3,3), Z(1:3));
% 
% locations = temp_data(1:4,:);
% locations = [locations; 1 0 0; 0 1 0; 0 0 1];
% A = locations(:, 1); B = locations(:, 2); C = locations(:, 3);
% Z = [1.0; 0.0; 0.0; 0.0;0;0;0]+0.1;
% ternpcolor(A,B,C,Z);





% locations = [1.000 0.000 0.000 ; 0.000 1.000 0.000; 0.000 0.000 1.000; 0.500 0.500 0.500];
% data = [1.0; 0.0; 0.0; 0.0]; A = locations(:, 1); B = locations(:, 2); C = locations(:, 3);
% terncontour(A, B, C, data)
% ternpcolor(A, B, C, data)