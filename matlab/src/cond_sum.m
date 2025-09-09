function est = cond_sum(normw,f,wf)
% computes approximate cancellation error based on condition number of
% summation operation (in the maximum norm)
kappa = normw*norm(f,inf)/abs(wf);
est = eps*kappa;
end
