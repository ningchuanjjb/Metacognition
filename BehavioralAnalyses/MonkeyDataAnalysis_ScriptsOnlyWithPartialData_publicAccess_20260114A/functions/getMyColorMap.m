function mycolormap = getMyColorMap()
mycolorpoint= [0   0   0;
               0   0   255;
               255 0   255; 
               255 0   0;
               255 255 0; 
               255 255 255]; 
cmin=0;
cmax=1;
clims=[cmin cmax]; % 定义颜色标尺范围
mycolorposition=linspace(clims(1),clims(2),6);
inp_100=linspace(cmin,cmax,100); % 插值数目100个
mycolormap_r=interp1(mycolorposition,mycolorpoint(:,1),inp_100,'linear','extrap');
mycolormap_g=interp1(mycolorposition,mycolorpoint(:,2),inp_100,'linear','extrap');
mycolormap_b=interp1(mycolorposition,mycolorpoint(:,3),inp_100,'linear','extrap');
mycolormap=[mycolormap_r',mycolormap_g',mycolormap_b']/256;
mycolormap=round(mycolormap*10^4)/10^4;%保留4位小数