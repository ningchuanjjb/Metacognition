% dftregistration_min_max
% dftregistration_min_max
% shift_reconstruct
% cell2mat_ov

close all
clear all
load temp_reg
% temp_imref2d = imref2d(size(M_final));
% [tform,peakcorr] = imregcorr(bw_template_2_shift,M_final);
% fprintf('peakcorr = %.2f.\n',peakcorr);
% I_reg = imwarp(bw_template_2_shift,tform,'OutputView',temp_imref2d);

PyramidLevels = 9;
iterNum = 1000;
temp_iterRange = zeros(1,PyramidLevels);
for tempi=1:PyramidLevels
    if tempi == 1
        temp_iterRange(tempi) = iterNum;
    else
        temp_iterRange(tempi) = ceil(temp_iterRange(tempi-1)/1.5);
    end
end
[D,~] = imregdemons(bw_template_2_shift,M_final,ceil(iterNum:-iterNum/PyramidLevels:1),...
    'AccumulatedFieldSmoothing',1,...
    'PyramidLevels',9,...%9
    'DisplayWaitbar',false);
I_reg = imwarp(bw_template_2_shift,D,'nearest');

% bw_template_2_shift_batch = true([size(bw_template_2_shift),2]);
% bw_template_2_shift_batch(:,:,1) = bw_template_2_shift;
% bw_template_2_shift_batch(:,:,2) = bw_template_2_shift;
% M_final_batch = true([size(M_final),2]);
% M_final_batch(:,:,1) = M_final;
% M_final_batch(:,:,2) = M_final;
% 
% numFrames = 2;
% % D_batch = repmat(D,[1,numFrames,1]);
% % D_batch = reshape(D_batch,[512,512,numFrames,2]);
% % D_batch(:,:,:,3) = 0;
% 
% % [D_batch,I_reg_batch] = imregdemons(bw_template_2_shift_batch,M_final_batch,ceil(iterNum:-iterNum/PyramidLevels:1),...
% %     'AccumulatedFieldSmoothing',0.5,...
% %     'PyramidLevels',9,...%9
% %     'DisplayWaitbar',false);
% 
% D_batch = zeros(512,512,numFrames,3);
% for tempi=1:3
%     if tempi < 3
%         for tempj=1:numFrames
%             D_batch(:,:,tempj,tempi) =  D(:,:,tempi);
%         end
%     elseif tempi == 3
%         %D_batch(:,:,:,tempi) = 0;        
%     end
% end
% 
% I_reg_batch = imwarp(bw_template_2_shift_batch,D_batch);
% % I_reg = mean(I_reg_batch,3);



I_1 = bw_template_1;
I_2 = M_final;
I_3 = bw_template_2_shift;
I_4 = I_reg;


fig1 = figure('Name','Fig1, I_1','NumberTitle','off');
set(gcf,'Position',[0 0 1920 1080]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(I_1);
colormap(gray);
axis off equal

fig2 = figure('Name','Fig2, I_2','NumberTitle','off');
set(gcf,'Position',[0 0 1920 1080]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(I_2);
colormap(gray);
axis off equal

fig3 = figure('Name','Fig3, I_3','NumberTitle','off');
set(gcf,'Position',[0 0 1920 1080]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(I_3);
colormap(gray);
axis off equal

fig4 = figure('Name','Fig4, I_4','NumberTitle','off');
set(gcf,'Position',[0 0 1920 1080]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(I_4);
colormap(gray);
axis off equal