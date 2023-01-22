function [Xhat min_reldiff best_alpha] = twg_rpca(X,W,Truth)
min_reldiff = Inf;
for alpha = 0.1:0.1:3
  result = wrapper_rpca(X,W,alpha);
  reldiff = norm(tensor(result.X)-tensor(Truth))/norm(tensor(Truth));
  if reldiff<min_reldiff
    Xhat = result.X;
    min_reldiff = reldiff;
    best_alpha = alpha;
  end
end
end