function X = cell2mat_ov_sum_jjb_v1(I,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap,sz,Bs)

% converts a cell array to a matrix when the cell elements overlap
% INPUTS:
% I:            cell array
% grid_size:    true size of each element
% overlap:      amount of overlap in each direction
% d1:           number of rows of matrix
% d2:           number of columns of matrix

% OUTPUT:
% X:            output matrix

% Written by Eftychios A. Pnevmatikakis, Simons Foundation, 2016

% t0 = tic;

if nargin < 10 || isempty(Bs)
    Bs = cellfun(@(x) ones(size(x)), I,'un',0);
end
X = zeros([sz,size(I{1,1},length(sz)+1)], 'single');
B = zeros(size(X), 'single');
if length(sz) == 2; sz(3) = 1; end


% toc(t0);
% a = 1;
% t1 = tic;

for i = 1:length(xx_f)
    for j = 1:length(yy_f)
        for k = 1:length(zz_f)            
            extended_grid = [max(xx_s(i)-overlap(1),1),min(xx_f(i)+overlap(1),sz(1)),max(yy_s(j)-overlap(2),1),min(yy_f(j)+overlap(2),sz(2)),max(zz_s(k)-overlap(3),1),min(zz_f(k)+overlap(3),sz(3))];
            
            I_ijk = I{i,j,k};
            Bs_ijk = Bs{i,j,k};
            
            Xtemp = zeros(size(I_ijk), 'single');
            Btemp = Xtemp;
            ind = ~isnan(I_ijk);            
            
            %t01 = tic;
            
            %Xtemp(ind) = Bs{i,j,k}(ind).*I{i,j,k}(ind);
            %Xtemp(ind) = Bs_ijk(ind).*I_ijk(ind);
            Xtemp = Bs_ijk.*I_ijk;
            Xtemp(~ind) = 0;
            
            %toc(t01);
            
            %Btemp(ind) = Bs_ijk(ind);
            Btemp = Bs_ijk;
            Btemp(~ind) = 0;
            
            %toc(t01);
            
            X(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) = ...
                X(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) + Xtemp; 
            
            %toc(t01);
            
            B(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) = ...
                B(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) + Btemp;            

            %toc(t01);
            %a = 1;
        end
    end
end
% toc(t1);
% a = 1;
% t2 = tic;

X = X./B;

% toc(t2);
% a = 1;

