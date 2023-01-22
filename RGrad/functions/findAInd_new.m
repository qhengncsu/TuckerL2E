function J = findAInd_new(G,alpha)
% input: G - the tensor
%        alpha - thersholding parameter
% output: J - a 0,1 tensor of size(G) and 1 if it is in the level-alpha 
%         active indice set, 0 otherwise

sz = size(G); 
m = numel(sz); 
n = prod(sz); 
dinv = n./sz;

J = ones(sz);
JJ = cell(m,1);
if alpha == 0 % exact low rank case without outliers
    return;
elseif alpha >= 1 
    J = ones(sz);
    return;
end

for i = 1:m
    k = ceil(dinv(i)*alpha);
    mat = ndim_unfold(G,i);
    JJ{i} = zeros(size(mat));
    
    for j = 1:sz(i)
        row = abs(mat(j,:));
        [~,I] = maxk(row,k);
        JJ{i}(j,I) = 1;
    end
    JJ{i} = ndim_fold(JJ{i},i,sz);
end

for i = 1:m
    J = J.*JJ{i};
end


