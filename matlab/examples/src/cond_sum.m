function est = cond_sum(normw,f,wf)
kappa = normw*norm(f,inf)/abs(wf);
est = eps*kappa;
end
