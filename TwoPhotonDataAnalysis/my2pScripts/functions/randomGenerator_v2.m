function randomValue = randomGenerator_v2(randNum, pdf, x_axis)

x_axis_step = x_axis(2) - x_axis(1);

cdf = pdf;
for tempi=1:length(cdf)
    cdf(tempi) = sum(pdf(1:tempi).*x_axis_step);
end


randomValue = zeros(1, randNum);

rand1 = rand(1,randNum);

rand1(rand1>max(cdf)) = max(cdf);

for tempi=1:randNum
    randomValue(tempi) = x_axis(find(cdf>=rand1(tempi), 1));
end

% a = 1;


