%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation of the ADMM in "Robust Tensor Decomposition with Gross Corruption"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,S] = tRPCA_new(A,lambda1,lambda2)
    if nargin < 2
        lambda1 = 1;
        lambda2 = 1.5e-2;
    end

    sz = size(A);
    m = length(sz);
    tol = 1e-6;
    maxit = 500;
    
    eta = 1e-4;
    rho = 1.1;
    max_eta = 1e4;
    
    psi = cell(m,1); % auxiliary variable
    a = cell(m,1); % Lagrange multiplier

    %% initialization
    T = zeros(sz);
    S = zeros(sz);
    for i = 1:m
        psi{i} = zeros(size(ndim_unfold(A,i)));
        a{i} = psi{i};
    end
    
    %% iteration
    for i = 1:maxit
        Tprev = T;
        % update T
        sum_a_psi = zeros(sz);
        for j = 1:m
            sum_a_psi = sum_a_psi + ndim_fold(a{j},j,sz) + eta * ndim_fold(psi{j},j,sz);
        end
        T = (-S + A + sum_a_psi)/(1+m*eta);
        
        % update S
        S = So(lambda2, A-T);
        
        % update psi, phi
        for j = 1:m
            psi{j} = Do(lambda1/(eta*m),ndim_unfold(T,j) - a{j}/eta);
        end
        
        for j = 1:m
            a{j} = a{j} + eta*(psi{j} - ndim_unfold(T,j));
        end
        
%         fprintf('Step: %d. Current error wrt previous step: %f, sparsity: %d\n',i,norm(T(:)-Tprev(:))/norm(T(:)),nnz(S));
        if i==1 || mod(i,10) == 0
            fprintf('Step: %d. Data fidelity: %f, r: (%d,%d,%d), l1 norm: %f, nuclear norm: (%f,%f,%f)\n',...
                i,norm(A(:)-T(:)-S(:)),rank(psi{1}),rank(psi{2}),rank(psi{3}),sum(abs(S(:))),sum(svd(psi{1})),sum(svd(psi{2})),sum(svd(psi{3})));
        end
            
        if norm(T(:)-Tprev(:))/norm(T(:)) < tol
            break;
        end
        eta = min(rho*eta,max_eta);
    end
end

function r = So(tau, X)
    % shrinkage operator
    r = sign(X) .* max(abs(X) - tau, 0);
end

function r = Do(tau, X)
    % shrinkage operator for singular values
    [U, S, V] = svd(X, 'econ');
    r = U*So(tau, S)*V';
end
        
        
        
