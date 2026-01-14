a = zeros(10,1024,'single','gpuArray')+1;
temp_cell = {a};
b = temp_cell{1};

temp_cell{2} = ones(10,1024,'single');

temp_sum = temp_cell{1}+temp_cell{2};



temp_struct = struct;
temp_struct.a = a;
c = temp_struct.a;
temp_struct.x = 2;
x = temp_struct.x;


