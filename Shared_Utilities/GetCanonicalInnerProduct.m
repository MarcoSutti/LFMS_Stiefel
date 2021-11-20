function [ inner_product ] = GetCanonicalInnerProduct( csi, zeta, Y0 )

% function [ inner_product ] = GetCanonicalInnerProduct( csi, zeta, Y0 )
% Purpose: Returns the canonical inner product of two tangent vectors 
%          belonging to the tangent space at X using the canonical metric
%          Stiefel manifold.

inner_product = trace( csi' * ( eye(size(Y0,1)) - 0.5 * (Y0 * Y0') ) * zeta );

end