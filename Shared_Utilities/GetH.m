function [ H ] = GetH( n, p )

% function H = GetH( n, p )
% Purpose: Returns the transformation matrix between blockwise
%          vectorization and ordinary vectorization.
% Created:     23.06.2016
% Last change: 30.06.2016

if n==p
    % if n=p, ordinary vectorization and blockwise vectorization coincide
    H = speye(n^2);
else
    % Build upper block
    for j=0:p-1
        H11( j*n+1:j*n+p, j*p+1:(j+1)*p ) = speye(p);
        H11(j*n+p+1:(j+1)*n, p^2+j*(n-p)+1:p^2+j*(n-p)+n-p ) = speye(n-p);
    end
    % Assemble into total H
    H( 1:n*p, 1:n*p ) = H11;
    
    % Build lower block
    for j=0:n-p-1
        H22( j*n+1:j*n+p, j*p+1:(j+1)*p ) = speye(p);
        H22( j*n+p+1:(j+1)*n, (n-p)*p+j*(n-p)+1:(n-p)*p+j*(n-p)+n-p ) = speye(n-p);
    end
    % Assemble into total H
    H( n*p+1:n^2, n*p+1:n^2 ) = H22;
end

end