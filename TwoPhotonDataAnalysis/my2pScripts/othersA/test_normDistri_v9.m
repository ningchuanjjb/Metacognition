clear
close all

% targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
% cd(targetPATH)

% x = normrnd(5,sqrt(4),1000,1);
% y = normrnd(5,sqrt(4),1000,1);
% z = x+y;
% mean(z), var(z)
% 
% mu = 10;
% sigma = sqrt(8);
% zmin = min(z); zmax = max(z); zrange = range(z);
% binw = zrange / 30;
% edges = zmin + binw*(0:30);
% n = histc(z, edges);
% bar(edges, n, 'histc');
% hold on
% xx = zmin:(zrange/1000):zmax;
% plot(xx, binw*length(z)*normpdf(xx,mu,sigma), 'r');
% title('Normal Fit');


a_range = [-2 3];
a_num = 5000;
x_mu = 0.3;
x_std = 0.05;
y_mu = 0.7;
y_std = 0.08;

a=linspace(a_range(1),a_range(2),a_num); 
a_step = a(2)-a(1);
b=linspace(2*a_range(1),2*a_range(2),a_num*2-1); 
x_pdf=normpdf(a,x_mu,x_std);
y_pdf=normpdf(a,y_mu,y_std); 
z_pdf=conv(x_pdf,y_pdf).*a_step;

weightX = 0.8;
weightY = 1 - weightX;% 2 - weightX;
% xy_pdf = (x_pdf.^weightX.*y_pdf.^weightY)./sum(x_pdf.^weightX.*y_pdf.^weightY.*a_step);
xy_pdf = (x_pdf.^weightX).*(y_pdf.^weightY)./sum((x_pdf.^weightX).*(y_pdf.^weightY).*a_step);

[~,I] = max(xy_pdf);
xy_mu = a(I);

xy_precise_std = x_std*y_std/sqrt(weightY*x_std^2+weightX*y_std^2);
temp_fun = @(t)fixed_point(t,x_mu,x_std,weightX,y_mu,y_std,weightY);
[~,mu3]=minsearch(temp_fun,mean([x_mu,y_mu]));
xy_precise_mu = (x_mu*y_std^2*weightX+y_mu*x_std^2*weightY)/(y_std^2*weightX+x_std^2*weightY);

xy_precise_pdf = normpdf(a,xy_precise_mu,xy_precise_std);


        %confidence_mu = weightX*x_mu + weightY*y_mu;
        %confidence_var = (weightX^2)*(x_std^2) + (weightY^2)*(y_std^2);
        %confidence_std = sqrt(confidence_var);
        

xlim_min = min( a(min(find(x_pdf>0.001,1), find(y_pdf>0.001,1))), b(find(z_pdf>0.001,1)) );
xlim_max = max( a(max(find(x_pdf>0.001,1,'last'), find(y_pdf>0.001,1,'last'))), b(find(z_pdf>0.001,1,'last')) );
ylim_min = 0;
% ylim_max = min(max([x_pdf,y_pdf,z_pdf])+0.05,1);
ylim_max = max([x_pdf,y_pdf,z_pdf,xy_pdf])+0.05;

figure(1)
plot(a,x_pdf,'Color',[0.75 0.25 0.25]);
hold on
plot(a,y_pdf,'Color',[0.25 0.75 0.25]);
hold on
plot(b,z_pdf,'Color',[0.25 0.25 0.75]);
hold on
plot(a,xy_pdf,'-','Color',[0.5 0.5 0.5],'LineWidth',5);
hold on
plot(a,xy_precise_pdf,'--','Color',[0 0 0]);
hold on
xlim([xlim_min,xlim_max]);
ylim([ylim_min,ylim_max]);


sum(x_pdf.*a_step);
sum(y_pdf.*a_step);
sum(z_pdf.*a_step);

x_cdf = x_pdf;
y_cdf = y_pdf;
for tempi=1:length(x_cdf)
    x_cdf(tempi) = sum(x_pdf(1:tempi).*a_step);
    y_cdf(tempi) = sum(y_pdf(1:tempi).*a_step);    
end
z_cdf = z_pdf;
for tempi=1:length(z_cdf)
    z_cdf(tempi) = sum(z_pdf(1:tempi).*a_step);
end
% figure(2)
% plot(a,x_cdf,'Color',[0.75 0.25 0.25]);
% hold on
% plot(a,y_cdf,'Color',[0.25 0.75 0.25]);
% hold on
% plot(b,z_cdf,'Color',[0.25 0.25 0.75]);
% hold on
% xlim([xlim_min,xlim_max]);

% randomGeneratorName = 'randomGenerator';
% randomGeneratorName_v = autoGetFunName_myScripts(randomGeneratorName, [targetPATH '\functions']);
% fprintf('Now runing %s.  ------> \n', randomGeneratorName_v);
% fun_randomGenerator = str2func(randomGeneratorName_v);
% 
% randNum = 50000;
% nbins = 50;
% FaceAlpha = 0.7;
% 
% x_random = fun_randomGenerator(randNum, x_pdf, a);
% y_random = fun_randomGenerator(randNum, y_pdf, a);
% z_random = fun_randomGenerator(randNum, z_pdf, b);
% xy_random = fun_randomGenerator(randNum, xy_pdf, a);
% 
% [muHat,sigmaHat] = normfit(xy_random);
% 
% temp_std = (x_std^weightX)*(y_std^weightY);

% figure(3)
% histogram(x_random,nbins,'Normalization','pdf','FaceColor',[0.75 0.25 0.25],'FaceAlpha',FaceAlpha);
% hold on
% histogram(y_random,nbins,'Normalization','pdf','FaceColor',[0.25 0.75 0.25],'FaceAlpha',FaceAlpha);
% hold on
% histogram(z_random,nbins,'Normalization','pdf','FaceColor',[0.25 0.25 0.75],'FaceAlpha',FaceAlpha);
% hold on
% xlim([xlim_min,xlim_max]);
% ylim([ylim_min,ylim_max]);



function  [out,t]=fixed_point(t,mu1,std1,w1,mu2,std2,w2)
out = std2^2*w1.*(t-mu1).^2 + std1^2*w2.*(t-mu2).^2;
a = 1; %#ok<*NASGU>
end

function [out,t]=minsearch(f,t_initial)
options = optimset('fminsearch');
options.TolFun = 1e-10;%1e-7-->1e-4-->1e-3
options.TolX = 1e-10;%1e-7-->1e-4-->1e-1
options.Display = 'off'; % off, iter, final

[t,out] = fminsearch(f, t_initial, options);
a = 1;
end

