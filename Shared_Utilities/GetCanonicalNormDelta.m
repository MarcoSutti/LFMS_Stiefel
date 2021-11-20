function [ norm_Delta ] = GetCanonicalNormDelta( Delta, Y0 )

% function [ norm_Delta ] = GetCanonicalNormDelta( Delta, Y0 )
% Purpose: Returns the canonical norm of Delta, which is defined on the
%          tangent space by using a quotient approach to the Stiefel manifold.

Omega = Y0'*Delta;

norm_Delta = sqrt( trace(Delta'*Delta) - 0.5*trace(Omega'*Omega) );

end