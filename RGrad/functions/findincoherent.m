function mu = findincoherent(U)

mu = 0;

l = length(U);
for i = 1:l
    V = U{i};
    sz = size(V);
    for j = 1:sz(1)
        m = norm(V(j,:))*sqrt(sz(1)/sz(2));
        if m >mu
            mu = m;
        end
    end
end