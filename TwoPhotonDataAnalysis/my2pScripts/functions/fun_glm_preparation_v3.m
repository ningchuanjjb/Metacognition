function glm_prep = fun_glm_preparation_v3(glm_prep_options)

order_glm = glm_prep_options.order_glm;
plot_lengthFlag = glm_prep_options.plot_lengthFlag;
sequence_lengthx_onehot = glm_prep_options.sequence_lengthx_onehot;
numFrames = glm_prep_options.numFrames;

order_glm_valid = min(order_glm,plot_lengthFlag);
sequence_lengthx_onehot_oneOrder = sequence_lengthx_onehot;
if order_glm_valid > 0
    for tempi=1:size(sequence_lengthx_onehot,1)
        a1 = find(sequence_lengthx_onehot(tempi,:)==true,order_glm_valid);
        a1 = a1(end);
        a2 = false(1,numFrames);
        a2(a1) = true;
        sequence_lengthx_onehot_oneOrder(tempi,:) = a2;
    end
end

sequence_lengthx_onehot_order = zeros(size(sequence_lengthx_onehot,1),plot_lengthFlag*numFrames);
for tempOrder=1:plot_lengthFlag
    temp_order_glm_valid = min(tempOrder,plot_lengthFlag);
    temp_onehot = sequence_lengthx_onehot;
    if temp_order_glm_valid > 0
        for tempi=1:size(sequence_lengthx_onehot,1)
            a1 = find(sequence_lengthx_onehot(tempi,:)==true,temp_order_glm_valid);
            a1 = a1(end);
            a2 = false(1,numFrames);
            a2(a1) = true;
            temp_onehot(tempi,:) = a2;
        end
    end
    temp_range = (1:numFrames) + numFrames*(tempOrder-1);
    sequence_lengthx_onehot_order(:,temp_range) = temp_onehot;
end

temp_locValidBool = true(1,numFrames);
if order_glm_valid > 0
    temp_locValidBool(numFrames-plot_lengthFlag+order_glm_valid+1:end) = false;
    temp_locValidBool(1:(order_glm_valid-1)) = false;
end

temp_locValidBool_real = true(1,numFrames);
for tempi=1:numFrames
    if sum(sequence_lengthx_onehot_oneOrder(:,tempi)) == 0
        temp_locValidBool_real(tempi) = false;
    end
end

glm_prep = struct;
glm_prep.order_glm_valid = order_glm_valid;
glm_prep.sequence_lengthx_onehot_oneOrder = sequence_lengthx_onehot_oneOrder;
glm_prep.sequence_lengthx_onehot_order = sequence_lengthx_onehot_order;
glm_prep.temp_locValidBool = temp_locValidBool;
glm_prep.temp_locValidBool_real = temp_locValidBool_real;