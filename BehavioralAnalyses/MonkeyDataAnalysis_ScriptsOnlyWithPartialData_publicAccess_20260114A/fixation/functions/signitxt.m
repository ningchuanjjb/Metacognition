function [signi] = signitxt(x1,y1,x2,y2,group1Num,group2Num,high)
a = [0.05,0.01,0.001];
signi = zeros(length(a),group1Num);
nn = group1Num * group2Num;
for i = 1:group1Num
    if length(y1{i}) <= 1 || length(y2{i}) <= 1
        continue
    else
        for ii = 1:length(a)
            alpha = a(ii)/nn;
            [x1_x2_sign,p2,ci2,stats2] = ttest2(y1{i},y2{i},'Alpha',alpha);
            signi(ii,i) = x1_x2_sign;
        end
    end
end
for i = 1:group1Num
    max_up = max([y1{i},y2{i}])+high+0.01;
    max_down = max([y1{i},y2{i}])+high;
    txt1_y = max_up + 0.01;
    if  isempty(find(signi(:,i), 1)) == 0 && sum(isnan(signi(:,i)))== 0
        line([x1(i),x2(i)],[max_up,max_up],'linestyle','-','color','k','linewidth',1);
        line([x1(i),x1(i)],[max_down,max_up],'linestyle','-','color','k','linewidth',1);
        line([x2(i),x2(i)],[max_down,max_up],'linestyle','-','color','k','linewidth',1);
    end
    if signi(3,i) == 1
        txt = '***';
        text((x1(i)-0.03+x2(i))/2-0.05,txt1_y,txt,'FontSize',18)
    elseif signi(2,i) == 1
        txt = '**';
        text((x1(i)-0.04+x2(i))/2-0.02,txt1_y,txt,'FontSize',18)
    elseif signi(1,i) == 1
        txt = '*';
        text((x1(i)-0.05+x2(i))/2,txt1_y,txt,'FontSize',18)
    end
end