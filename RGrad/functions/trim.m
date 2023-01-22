function trimU = trim(U,mu0)

sz = size(U);
trimU = zeros(sz);
for i = 1:sz(1)
    if norm(U(i,:)) > sqrt(mu0*sz(2)/sz(1))
        trimU(i,:) = U(i,:)*sqrt(mu0*sz(2)/sz(1))/norm(U(i,:));
    else
        trimU(i,:) = U(i,:);
    end
end
        