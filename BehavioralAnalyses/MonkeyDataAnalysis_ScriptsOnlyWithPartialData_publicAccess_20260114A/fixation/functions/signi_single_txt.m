function [signi] = signi_single_txt(x1,y1,X,group1Num,group2Num,high)
type = isa(y1,'double');
if type == 0
    a = [0.05,0.01,0.001];
    signi = zeros(length(a),group1Num);
    nn = group1Num * group2Num;
    for i = 1:group1Num
        if length(y1{i}) <= 1
            continue
        else
            for ii = 1:length(a)
                alpha = a(ii)/nn;
                [x1_x2_sign,p2,ci2,stats2] = ttest(y1{i},X,'Alpha',alpha);
                signi(ii,i) = x1_x2_sign;
            end
        end
    end
    for i = 1:group1Num
        max_up = max(y1{i})+high+0.01;
        txt1_y = max_up + 0.01;
        if signi(3,i) == 1
            txt = '***';
            text(x1(i)-0.1,txt1_y,txt,'FontSize',18)
        elseif signi(2,i) == 1
            txt = '**';
            text(x1(i)-0.05,txt1_y,txt,'FontSize',18)
        elseif signi(1,i) == 1
            txt = '*';
            text(x1(i),txt1_y,txt,'FontSize',18)
        end
    end
    
elseif type == 1
    a = [0.05,0.01,0.001];
    signi = zeros(length(a),1);
    nn = group1Num * group2Num;
    for ii = 1:length(a)
        alpha = a(ii)/nn;
        [x1_x2_sign,p2,ci2,stats2] = ttest(y1,X,'Alpha',alpha);
        signi(ii) = x1_x2_sign;
    end
    max_up = max(y1)+high+0.01;
    txt1_y = max_up + 0.01;
    if signi(3) == 1
        txt = '***';
        text(x1-0.1,txt1_y,txt,'FontSize',18)
    elseif signi(2) == 1
        txt = '**';
        text(x1-0.05,txt1_y,txt,'FontSize',18)
    elseif signi(1) == 1
        txt = '*';
        text(x1,txt1_y,txt,'FontSize',18)
    end
end
