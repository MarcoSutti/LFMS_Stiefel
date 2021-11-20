function [ J_Z1_x_FD ] = GetJZ1xFD( E, Q, x, h )

% function [ J_Z1_x_FD ] = GetJZ1xFD( E, Q, x, h )
% Purpose: Returns the Jacobian J_Z1_x used in the single shooting method
%          on the Stiefel manifold, calculated using Finite Differences.
%          We use the formulation which expresses Z1 as a function of x.
%          x is the column vector containing the free parameters that build
%          up a tangent vector in the tangent space T_{Y0}St(n,p).
%          See notes of 21.09.2016.
%          See also the scripts Derivative_Z1x_wrt_x and
%          Agreement_J_Z1_x_analytic_and_J_Z1_x_FD, in which we
%          numerically verified the correctness of these formulas.
% Created:     22.09.2016
% Last change: 22.09.2016

n = size(Q,1);
p = sqrt(size(E,1));

length_x = length(x);

%--------------------------------------------------------------------------
% Initialize matrices that depend only on n and p
I_length_x = eye(length_x);
%--------------------------------------------------------------------------

% Build Jacobian with finite differences:
for k=1:length_x
    e_k = I_length_x(:,k);
    Z1_x_plus_hek  = GetZ1x( E, Q, x + h * e_k );
    Z1_x_minus_hek = GetZ1x( E, Q, x - h * e_k );
    dZ1dx = ( Z1_x_plus_hek - Z1_x_minus_hek )/(2*h);
    J_Z1_x_FD(:,k) = dZ1dx(:);   
end

end