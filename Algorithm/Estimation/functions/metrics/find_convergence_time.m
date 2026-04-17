function tc = find_convergence_time(t, err, threshold)
idx = find(err <= threshold, 1, 'first');
if isempty(idx)
    tc = inf;
else
    tc = t(idx);
end
end
