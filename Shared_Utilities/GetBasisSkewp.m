function [ E ] = GetBasisSkewp( p )

% function [ E ] = GetBasisSkewp( p )
% Purpose: Returns E, a basis in vectorized form for space of p-by-p
%          skew-symmetric matrices.
%          This basis E allows to get a vectorization of a skew-symmetric
%          matrix Omega as vecOmega = E*s, where s is a column vector of
%          (1/2)*p*(p-1) random parameters.
% Created:     22.09.2016
% Last change: 31.10.2016

if p==1
    E = 1; % If p=1, we do not have a skew-symmetric matrix, but a scalar, equal to zero.
else
    Ip = eye(p);
    
    % Build basis for Skew(p):
    E = [];
    
    for i=2:p
        for j=1:i-1
            E = [ E, reshape( Ip(:,i)*Ip(j,:)-Ip(:,j)*Ip(i,:), [p^2,1] ) ];
        end
    end
    
    % NB: We can also enumerate the bases Eij along the columns. In this case,
    % change the loop with for j=1:p-1
    %                          for i=j+1:p
end
end
