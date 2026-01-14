temptemp_path = 'D:\Downloads\20231212A';
temptemp_loadPic = 'pic6.png';
temptemp_sacePic = 'pic6_gray.png';

I_load = imread([temptemp_path,'\',temptemp_loadPic]);

image(I_load);

I_save = rgb2gray(I_load);

image(I_save);

imwrite(I_save,[temptemp_path,'\',temptemp_sacePic]);