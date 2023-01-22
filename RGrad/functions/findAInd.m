function J = findAInd(G,alpha)
% input: G - the tensor
%        alpha - thersholding parameter
% output: J - a 0,1 tensor of size(G) and 1 if it is in the level-alpha 
%         active indice set, 0 otherwise

sz = size(G); 
m = numel(sz); 
n = prod(sz); 
dinv = n./sz;
thres = cell(m,1);
tol = cell(m,1);

J = zeros(sz);
if alpha == 0 % exact low rank case without outliers
    return;
elseif alpha >= 1 
    J = ones(sz);
    return;
end
% pre-compute the threshold
tick = zeros(sz);
for i = 1:m
    thres{i} = zeros(sz(i),1);
    k = ceil(dinv(i)*alpha);
    tol{i} = k*ones(sz(i),1);
    mat = ndim_unfold(G,i);
    unfold_tick = ndim_unfold(tick,i);
    
%     absmat = abs(mat);
%     [kmax,~] = maxk(absmat,k,2);
%     thres{i} = kmax(:,end);
%     for j = 1:sz(i)
%         unfold_tick(j,absmat(j,:)>=thres{i}(j)) = 1;
%     end
    
    for j = 1:sz(i)
        row = abs(mat(j,:));
        [kmax,~] = maxk(row,k);
        thres{i}(j) = kmax(end);
        unfold_tick(j,row>=thres{i}(j)) = 1;
    end
    tick = ndim_fold(unfold_tick,i,sz);
end


ind = cell(m,1);
for i = 1:n
    if tick(i) == 0
        continue;
    end
    
    sign = 0;
    [ind{:}] = ind2sub(sz,i);
    for k = 1:m
        if tol{k}(ind{k}) <= 0
            break;
        end
        
        if abs(G(ind{:})) < thres{k}(ind{k})
            break;
        end
        
        if k == m
            sign = 1;
        end
    end
    
    if sign == 1 
        J(ind{:}) = 1;
        for k = 1:m
            tol{k}(ind{k}) = tol{k}(ind{k}) - 1;
        end
    end
end

