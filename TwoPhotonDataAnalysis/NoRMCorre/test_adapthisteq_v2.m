clear all
% close all

existingTemplate_fileName = 'template_20230122A_Ding_Site09B_20230202T101726.tif';
existingTemplate = single(read_file(existingTemplate_fileName));
max_template = max(existingTemplate,[],'all');
I = existingTemplate./max_template;
% I = existingTemplate;
J1 = I.^0.4;
% J2 = adapthisteq(I,'Range','original');
J2 = adapthisteq(I,'NBins',512,'Range','original');
% J3 = adapthisteq(I,'NBins',1024,'Range','original','Distribution','rayleigh');
J3 = adapthisteq(I,'NBins',4096,'Range','original','Distribution','rayleigh','Alpha',0.6);
% J3 = adapthisteq(I,'NBins',4096,'Distribution','rayleigh','Alpha',0.6);
J4 = adapthisteq(I,'NBins',4096,'Range','original','Distribution','rayleigh','Alpha',0.8);
% J4 = adapthisteq(I,'NBins',4096,'Range','original','Distribution','rayleigh','Alpha',0.6,'ClipLimit',1);

% J4 = adapthisteq(I,'NBins',4096,'Distribution','rayleigh','Alpha',0.6);
% for tempi=1:10
%    %J4 = adapthisteq(J4,'NBins',4096,'Distribution','rayleigh','Alpha',0.6); 
%    J4 = adapthisteq(I,'NBins',4096,'Range','original','Distribution','rayleigh','Alpha',0.6,'ClipLimit',1);
% 
% end

% B = rescale(A,l,u)

fig1 = figure(1);
set (gca,'position',[0.01,0.01,0.98,0.98])
% imagesc(I.^(0.4));
% imagesc(J1);
imshow(J1);
colormap(gray);
axis off equal
fig12 = figure(12);
imhist(J1);


fig2 = figure(2);
set (gca,'position',[0.01,0.01,0.98,0.98])
% imagesc(J2);
imshow(J2);
colormap(gray);
axis off equal
fig22 = figure(22);
imhist(J2);

fig3 = figure(3);
set (gca,'position',[0.01,0.01,0.98,0.98])
% imagesc(J3);
imshow(J3);
colormap(gray);
axis off equal
fig32 = figure(32);
imhist(J3);

fig4 = figure(4);
set (gca,'position',[0.01,0.01,0.98,0.98])
% imagesc(J4);
imshow(J4);
colormap(gray);
axis off equal
fig42 = figure(42);
imhist(J4);
