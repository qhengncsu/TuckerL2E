function [T] = lowrankrecover(A)

    sz = size(A);
    m = length(sz);
    tol = 1e-10;
    maxit = 200;
    
    eta = 1e-1;
    rho = 1.1;
    max_eta = 1e4;
    
    lambda1 = 0.001;
    
    psi = cell(m,1); % auxiliary variable
    a = cell(m,1); % Lagrange multiplier

    %% initialization
    T = zeros(sz);
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
        T = (A + sum_a_psi)/(1+m*eta);
        
        % update psi
        for j = 1:m
            psi{j} = Do(lambda1/(eta),ndim_unfold(T,j) - a{j}/eta);
        end
        
        for j = 1:m
            a{j} = a{j} + eta*(psi{j} - ndim_unfold(T,j));
        end
        
%         fprintf('Step: %d. Current error wrt previous step: %f, sparsity: %d\n',i,norm(T(:)-Tprev(:))/norm(T(:)),nnz(S));
        if i==1 || mod(i,10) == 0
            fprintf('Step: %d. Data fidelity: %f, r: (%d,%d,%d), nuclear norm: (%f,%f,%f)\n',...
                i,norm(A(:)-T(:)),rank(psi{1}),rank(psi{2}),rank(psi{3}),sum(svd(psi{1})),sum(svd(psi{2})),sum(svd(psi{3})));
        end
            
        if norm(T(:)-Tprev(:))/norm(T(:)) < tol
            fprintf('Step: %d. Data fidelity: %f, r: (%d,%d,%d), nuclear norm: (%f,%f,%f)\n',...
                i,norm(A(:)-T(:)),rank(psi{1}),rank(psi{2}),rank(psi{3}),sum(svd(psi{1})),sum(svd(psi{2})),sum(svd(psi{3})));
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
        
        
        
