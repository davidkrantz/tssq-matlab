function [u1,u2,u3,stats] = ssq_sbt(curve, density, targets, opts)
%SSQ_SBT  Standard Singularity Swap Quadrature (SSQ) for slender-body 
% theory (SBT) kernel.
%
% Calls the same engine as TSSQ but with opts.use_mod = false.

    opts.use_mod = false;
    [u1,u2,u3,stats] = tssq_sbt(curve,density,targets,opts);
end
