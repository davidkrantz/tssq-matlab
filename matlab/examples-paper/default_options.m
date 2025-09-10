function opts = default_options(basis)
%DEFAULT_OPTIONS  Default options for SSQ/TSSQ bases used in examples.
%
%   opts = DEFAULT_OPTIONS(basis)
%   Returns a struct with default settings for the basis choice.
%
%   Different fields:
%     tol                 tolerance
%     nquad               quadrature order per panel 
%     rho                 Berstein radius limit rule for special quadrature
%     upsample            whether to upsample Gauss nodes
%     slender_eps         slender body radius
%     use_bjorck_pereyra  use Björck–Pereyra algorithm for Vandermonde solves
%     use_mod             use modified basis (TSSQ) if true, standard SSQ if false
%     corrR3              boolean, corect 1/R^3 with TSSQ
%     corrR5              boolean, correct 1/R^5 with TSSQ
%     basis               'monomial' = panel Gauss–Legendre (open curves) or 'fourier'  = global trapezoidal (closed curves)
%     Hlim                sets how much the adaptive quadrature refines

switch lower(basis)
    case 'monomial'
        opts = struct();
        opts.tol                = 1e-4;
        opts.nquad              = 16;
        opts.rho                = 3;
        opts.upsample           = true;
        opts.slender_eps        = 1e-3;
        opts.use_bjorck_pereyra = true;
        opts.use_mod            = true;
        opts.corrR3             = true;
        opts.corrR5             = true;
        opts.basis              = 'monomial';
        opts.Hlim               = 1;
    case 'fourier'
        opts = struct();
        opts.tol                = 1e-10;
        opts.nquad              = 16; % for adaptive quadrature reference
        opts.slender_eps        = 1e-3;
        opts.use_mod            = true;
        opts.corrR3             = true;
        opts.corrR5             = true;
        opts.basis              = 'fourier';
        opts.Hlim               = 4; % for adaptive quadrature reference
end
