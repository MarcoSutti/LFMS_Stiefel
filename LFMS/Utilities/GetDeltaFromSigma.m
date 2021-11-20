function [ length_Delta, Delta ] = GetDeltaFromSigma( n, p, m, Sigma_LF_k )

% function [ length_Delta, Delta ] = GetDeltaFromSigma( n, p, m, Sigma_LF_k )
% Purpose: Recover the tangent vector at Y0 to St(n,p) from the Sigma.
% Created: 15.02.2017.
% Last change: 15.02.2017.

% Find norms of tangent vectors from Sigma12 to Sigma(m-1,2)
Norm_Deltas = zeros(1,m-1);
for k=1:m-1
    Norm_Deltas(k) = GetCanonicalNormDelta( GetSigmakj( Sigma_LF_k, k, 2, n, p, m ), GetSigmakj( Sigma_LF_k, k, 1, n, p, m ) );
end

% These lengths should sum up to the norm of the original tangent vector, Delta_exact
length_Delta = sum(Norm_Deltas);

Delta = length_Delta * ( GetSigmakj( Sigma_LF_k, 1, 2, n, p, m )/Norm_Deltas(1) );

end