function k = get_k_vec(M,L)
% GET_K_VEC returns the wavenumbers associated with the equidistant grid of
%   length L discretized using M points. Takes into account if M is even or
%   odd.
% 
% k = get_k_vec(M,L) returns the wavenumbers in a row vector format
%
% INPUTS:
%   M - number of points on grid
%   L - lenght of interval, typically 2*pi
%
% OUTPUT:
%   k - row vector containing the wavenumbers 
% 
% AUTHOR: https://github.com/ludvigak/linequad
if mod(M,2)==0
    MM = M/2;
    k = (2*pi/L)*(-MM:(MM-1));
elseif mod(M-1,2)==0
    MM = (M-1)/2;
    k = (2*pi/L)*(-MM:MM);
else 
    error("k-vectors not computed");
end
end
