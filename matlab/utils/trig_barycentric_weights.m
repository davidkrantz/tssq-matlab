function w = trig_barycentric_weights(t, a)
%TRIG_BARYCENTRIC_WEIGHTS  Lagrange weights for periodic trigonometric interpolation
%
%  Given n equispaced nodes t in [0, 2*pi) and an evaluation point a,
%  returns column vector w0 such that sigma(a) ≈ sigma.' * w0.
%
%  Works for both even and odd n.

n = numel(t);
t = t(:);

% wrap a and differences to principal interval
a = mod(a, 2*pi);
d = mod(a - t + pi, 2*pi) - pi;   % differences in (-pi, pi]

% exact hit: Kronecker delta
tol = 1e-14;
hit = find(abs(d) < tol, 1);
if ~isempty(hit)
    w = zeros(n,1);
    w(hit) = 1;
    return
end

half = d/2;
nh   = (n/2)*d;

if mod(n,2)==1
    % n odd: l_j(a) = (1/n) * sin(n d / 2) / sin(d/2)
    w = (1/n) * ( sin(nh) ./ sin(half) );
else
    % n even: l_j(a) = (1/n) * sin(n d / 2) * cot(d/2)
    %        = (1/n) * sin(n d / 2) .* (cos(d/2) ./ sin(d/2))
    w = (1/n) * ( sin(nh) .* (cos(half) ./ sin(half)) );
end

w = w(:);
w = w / sum(w);
end
