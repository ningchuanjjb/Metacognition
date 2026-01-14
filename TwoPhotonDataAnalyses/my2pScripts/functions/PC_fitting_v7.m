function y = PC_fitting_v7(coeff,x,temp_boolIndex_location,fit_PC_num,fit_linear1_quadratic2)

a = 1;

numFrames = size(temp_boolIndex_location,1);

% locDis = zeros(size(x)); %#ok<*NASGU>
locDis = zeros(numFrames,1); %#ok<*NASGU>
for tempi=1:numFrames
    temp_disSquare = (x(tempi,:)-x).^2;
    
    temp_disSquare(tempi,:) = [];
    
    temp_disSquare2 = sum(temp_disSquare,2);
    temp_dis = sqrt(min(temp_disSquare2));
    
    %temp_dis = sqrt(sum(temp_disSquare,1));
    %temp_dis = sqrt(sum(temp_disSquare,'all'));
    locDis(tempi,:) = temp_dis;    
end
locDis_n11n = locDis/sum(locDis);

% locDis_n11n = ones(numFrames,1);

temp_locProb = zeros(numFrames,1);
for tempi=1:numFrames
    if fit_PC_num == 1
        if fit_linear1_quadratic2 == 1
            temp_locProb(tempi) = ...
                coeff(1)*x(tempi,1) ...
                +coeff(2);
        elseif fit_linear1_quadratic2 == 2
            temp_locProb(tempi) = ...
                coeff(1)*x(tempi,1).^2+coeff(2)*x(tempi,1) ...
                +coeff(3);
        end
    elseif fit_PC_num == 2
        if fit_linear1_quadratic2 == 1
            temp_locProb(tempi) = ...
                coeff(1)*x(tempi,1) ...
                +coeff(2) ...
                +coeff(3)*x(tempi,2);               
        elseif fit_linear1_quadratic2 == 2
            temp_locProb(tempi) = ...
                coeff(1)*x(tempi,1).^2+coeff(2)*x(tempi,1) ...
                +coeff(3) ...                
                +coeff(4)*x(tempi,2).^2+coeff(5)*x(tempi,2);
        end
    end
end
% temp_locProb(temp_locProb<0) = 0;
% temp_locProb(temp_locProb>1) = 1;

a = 1;


seqNum = size(temp_boolIndex_location,2);

temp_seqProb = zeros(seqNum,1);
for tempi=1:seqNum
    tempSeq_bool = temp_boolIndex_location(:,tempi);
    %temp_seqProb(tempi) = prod(temp_locProb(tempSeq_bool)) * prod(1-temp_locProb(~tempSeq_bool));
    temp_seqProb(tempi) = prod(temp_locProb(tempSeq_bool));
end

% temp_seqProb(temp_seqProb<0) = 0;
% temp_seqProb(temp_seqProb>1) = 1;

y = temp_seqProb;

a = 1;
