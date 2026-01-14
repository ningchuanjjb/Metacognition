clear
for i = 1:5
    %f(i) = parfeval(backgroundPool,@rand, 1, 1);
    F_A(i) = parfeval(backgroundPool,@funA, 2, 3,4,5);
    F_B(i) = afterEach(F_A(i),@funB,2);
end
% afterEach(f,@disp,0);

% F_B = afterEach(F_A,@funB,2);
% [A_output1,A_output2] = fetchOutputs(F_A);
% [B_output1,B_output2] = fetchOutputs(F_B);


function [out1,out2]=funA(in1,in2,in3)

out1=in1+in2+in3;
out2=in1-in2+in3;
end

function [out1,out2]=funB(in1,in2)

out1=max(in1,in2)^2;
out2=min(in1,in2)^2;
end