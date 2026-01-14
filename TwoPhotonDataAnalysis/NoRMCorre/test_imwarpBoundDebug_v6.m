
col_shift_multi = 10*ones(1,lY,'gpuArray');
row_shift_multi = 0*ones(1,lY,'gpuArray');

col_InfShift = 512*lY*2;

displacementField_1d_col = col_shift_multi * -1;
displacementField_1d_col_int = ceil(abs(displacementField_1d_col)).*sign(displacementField_1d_col);
col_leftBound_multi = max(1,1 - displacementField_1d_col_int);
col_rightBound_multi = min(512,512-displacementField_1d_col_int);

temp_grid_col = 1:512;
temp_grid_col = repmat(temp_grid_col,[512 1]);
temp_grid_col2 = reshape(temp_grid_col',[1 1 512 512]);
col_leftBound_multi2 = col_leftBound_multi';
col_rightBound_multi2 = col_rightBound_multi';
temp_bool = temp_grid_col2>=col_leftBound_multi2 & temp_grid_col2<=col_rightBound_multi2;

displacementField_3d_col = repmat(col_shift_multi' * -1,[1,1,512,512]);
displacementField_3d_col_valid = displacementField_3d_col.*temp_bool + col_InfShift.*(~temp_bool);

displacementField_3d_row = repmat(row_shift_multi' * -1,[1,1,512,512]);

displacementField = [displacementField_3d_col_valid,displacementField_3d_row];
displacementField = permute(displacementField, [4 3 1 2]);
displacementField = reshape(displacementField,[512,512*lY,2]);



% Shift reconstruct
Mf = reshape(Ytm,[512,512*lY]);
%Mf = imwarp(Mf,gather(displacementField),'linear');
Mf = imwarp(Mf,displacementField,'linear');
Mf = reshape(Mf,512,512,lY);

temp_Ytm_mean = gather(mean(Ytm,3));
temp_Ytm_mean = temp_Ytm_mean./max(temp_Ytm_mean,[],'all');
temp_Mf_mean = gather(mean(Mf,3));
temp_Mf_mean = temp_Mf_mean./max(temp_Mf_mean,[],'all');

close all

figure(1);
set(gcf,'Position',[35+400 35+0 850 850]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98])
imshow(temp_Ytm_mean.^0.4);
colormap(gray);
axis off equal

figure(2);
set(gcf,'Position',[35+400 35+0 850 850]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98])
imshow(temp_Mf_mean.^0.4);
colormap(gray);
axis off equal




