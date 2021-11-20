function [ Q, R ] = qrPosDiagR( X )

% function [ Q, R ] = qrPosDiagR( X )
% Purposes: Forces the diagonal of the R factor of the thin qr factorization to
%           be positive.

[ Q1, R1 ] = qr( X, 0 );   % get the orthogonal factor of X

% Force the diagonal of R to be positive:
R = diag( sign(diag(R1)) )*R1;
Q = Q1*diag( sign(diag(R1)) );

end
