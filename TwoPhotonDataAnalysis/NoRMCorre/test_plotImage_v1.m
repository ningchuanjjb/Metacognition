
n = 0;

n = n + 1;
figure(n);
set(gcf,'Position',[35+400 35+0 850 850]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98])
J = template_in;
J = J ./ max(template_in,[],'all');
imshow(J.^0.4);
colormap(gray);
axis off equal