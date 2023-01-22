function [T,S] = tRPCA(A)

    sz = size(A);
    m = length(sz);
    dprod = numel(A);
    tol = 1e-5;
    maxit = 500;
    
    eta = 1;
    rho = 1.1;
    max_eta = 1e4;
    
    lambda1 = 0.1;
    lambda2 = 0.01;
    
    psi = cell(m,1); % auxiliary variable
    phi = cell(m,1); % auxiliary variable
    a = cell(m,1); % Lagrange multiplier
    b = cell(m,1); % Lagrange multiplier

    %% initialization
    T = A;
    S = zeros(sz);
    for i = 1:m
        phi{i} = zeros(size(ndim_unfold(A,i)));
        psi{i} = phi{i};
        a{i} = phi{i};
        b{i} = a{i};
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
        sum_b_phi = zeros(sz);
        for j = 1:m
            sum_b_phi = sum_b_phi + ndim_fold(b{j},j,sz) + eta * ndim_fold(phi{j},j,sz);
        end
        S = (-T + A + sum_b_phi)/(1+m*eta);
        
        % update psi, phi
        for j = 1:m
            psi{j} = Do(lambda1/(eta*m),ndim_unfold(T,j) - a{j}/eta);
        end
        for j = 1:m
            phi{j} = So(lambda2/(eta*m),ndim_unfold(S,j) - b{j}/eta);
        end
        
        for j = 1:m
            a{j} = a{j} + (psi{j} - ndim_unfold(T,j));
            b{j} = b{j} + (phi{j} - ndim_unfold(S,j));
        end
        
        fprintf('Step: %d. Current error wrt previous step: %f\n',i,norm(T(:)-Tprev(:))/norm(T(:)));
        
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
        
        
        
