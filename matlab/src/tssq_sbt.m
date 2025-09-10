function [u1,u2,u3,stats] = tssq_sbt(curve, density, targets, opts)
%TSSQ_SBT  Translated Singularity Swap Quadrature (TSSQ) for slender-body 
% theory (SBT) kernel.
%
%   [u1,u2,u3,stats] = TSSQ_SBT(curve, density, targets, opts)
%
%   High-level user-facing driver that evaluates the Stokes slender-body
%   layer potential close to a filament using either panel-based 
%   Gauss-Legendre quadrature (monomial engine) or the global trapezoidal
%   rule (Fourier engine).
%
%   STOKES SLENDER-BODY KERNEL:
%       The potential at a target point x is
%
%           u(x) = ∫_Γ [  S(x - y) + (varrho^2 / 2) D(x - y) ] sigma(y) ds(y),
%
%       where sigma(y) is the force density along the filament, r = x - y, and
%       the kernels are the classical Stokeslet and stresslet doublet:
%
%           S(r) = I/|r| + (r r^T)/|r|^3,
%           D(r) = I/|r|^3 - 3 (r r^T)/|r|^5.
%
%       Here varrho = opts.slender_eps is the filament radius (regularization).
%
%       This splitting produces three groups of terms (1/R, 1/R^3, 1/R^5)
%       which are each evaluated using SSQ or TSSQ.
%
%   INPUT:
%     curve   : struct of function handles defining the geometry
%                 .x(t), .y(t), .z(t)  : position (R^3) for t in [0,1] or [0,2π)
%                 .s(t)                : speed |γ'(t)|
%               (optionally also .xp, .yp, .zp for user scripts)
%
%     density : struct of function handles for the line density
%                 .f1(t), .f2(t), .f3(t)
%
%     targets : 3 x Nt array of target points [X; Y; Z]
%
%     opts    : options struct. Defaults are supplied by DEFAULT_OPTIONS. 
%               Relevant fields are:
%                 .basis   : 'monomial' or 'fourier'
%                 .tol     : adaptive discretization tolerance (monomial path)
%                 .nquad   : #GL nodes per panel (monomial path)
%                 .rho     : limit-rule parameter for near weights (monomial)
%                 .upsample: logical, whether to upsample panel nodes (monomial)
%                 .slender_eps : SBT regularization parameter
%                 .use_bjorck_pereyra : use Björck–Pereyra Vandermonde solve
%                 .use_mod : true => TSSQ (modified/translated basis);
%                            false => SSQ (standard)
%                 .corrR3, .corrR5 : (Fourier/monomial engines) whether to
%                                    apply R^{-3} and/or R^{-5} corrections
%                                    in the translated basis evaluation
%
%   OUTPUT:
%     u1,u2,u3 : 1 x Nt component arrays of the velocity/potential
%     stats    : struct of timings and counters (engine-dependent fields)
%
%   NOTES:
%   * 'monomial' path:
%       - Builds an adaptive Gauss–Legendre panelization using ADAPTIVE_PANELIZATION
%         and calls INTERPOLATORY_QUADRATURE_MONOMIAL, which implements
%         SSQ/TSSQ in the monomial/translated monomial basis.
%   * 'fourier' path:
%       - Builds a global trapezoidal grid (via ADAPTIVE_GLOBAL_DISCRETIZATION
%         used here with tight tolerance to produce an equispaced grid) and calls
%         INTERPOLATORY_QUADRATURE_FOURIER, which implements SSQ/TSSQ in the
%         Fourier/modified Fourier basis.
%
%   AUTHOR: David Krantz (davkra@kth.se)

    if nargin < 4 || isempty(opts)
        opts = default_options_filament();
    end

    switch lower(opts.basis)
        case 'monomial'
            [tj, wj, npan, edges] = adaptive_panelization(curve.s, opts.nquad, opts.tol);
            fprintf('tol=%.1e, nquad=%d, npan=%d\n', opts.tol, opts.nquad, npan)
        
            xj = curve.x(tj); yj = curve.y(tj); zj = curve.z(tj);
            sj = curve.s(tj);
            f1j = density.f1(tj); f2j = density.f2(tj); f3j = density.f3(tj);
        
            [u1,u2,u3,~,~,stats] = interpolatory_quadrature_monomial( ...
                xj,yj,zj,sj,tj,wj,f1j,f2j,f3j, ...
                targets(1,:),targets(2,:),targets(3,:), ...
                curve.x,curve.y,curve.z,curve.s, ...
                opts.nquad,edges,opts.rho,opts.upsample, ...
                opts.slender_eps,opts.tol,opts.use_bjorck_pereyra,opts.use_mod, ...
                opts.corrR3, opts.corrR5);
        case 'fourier'
            [tj,wj] = adaptive_global_discretization(curve.s,5e-14);
            nquad_global = numel(tj);
            fprintf('tol=%.1e, nquad_global=%d\n', opts.tol, nquad_global);

            xj = curve.x(tj); yj = curve.y(tj); zj = curve.z(tj);
            sj = curve.s(tj);
            f1j = density.f1(tj); f2j = density.f2(tj); f3j = density.f3(tj);
            
            [u1,u2,u3,~,~,stats] = interpolatory_quadrature_fourier(...
                xj,yj,zj,sj,tj,wj,f1j,f2j,f3j, ...
                targets(1,:),targets(2,:),targets(3,:), ...
                curve.x,curve.y,curve.z,curve.s, ...
                nquad_global,opts.slender_eps,opts.tol,opts.use_mod, ...
                opts.corrR3, opts.corrR5);
    end
end
