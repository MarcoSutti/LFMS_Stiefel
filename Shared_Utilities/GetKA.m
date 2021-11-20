function [ KA ] = GetKA( A )

% function [ KA ] = GetKA( A )
% Purpose: Returns the Kronecker representation of the Fréchet derivative
%          of the matrix exponential at A.
%          From "Functions of Matrices", Higham, p.238, eq. (10.17b).
%          Originally from "Derivatives of the Matrix Exponential and Their
%          Computation", 1995, Najfeld and Havel, eq. (109).
% Created:     01.07.2016
% Last change: 31.08.2016

n = size(A,1);

B = (1/2)*( kron( A', eye(n) ) - kron( eye(n), A ) );

%lambda_max_A = max(eig(A))

Sinch = MatSinch( B );
% eig(B)
% eig(Sinch)
% eig(A)

% The differential of the matrix exponential in Kronecker representation (10.17b):
KA = kron( expm(A'/2), expm(A/2) ) * Sinch;

end


function [ Sinch ] = MatSinch( B )

% function [ Sinch ] = MatSinch( B )
% Purpose: Returns the matrix sinch function of B.
% Created:     01.07.2016
% Last change: 20.01.2017

% We define Sinch as a matrix function Sinch = V*diag(sinh(diag(D))./diag(D))/V;
% TODO: A skew-symmetric matrix should be diagonalised through a unitary
% matrix. MATLAB eig does not compute this. We should address this issue
% sooner or later.
[ V, D ] = eig(B);

d = imag(diag(D));

Sinch = V*diag(sinc(d/pi))/V;

% 30.08.2016:
% Determine whether B is diagonalizable. If B is nondiagonalizable, we
% have to use another formula.
% if rank(V)~= n2
%     error('B is not diagonalizable!!!') %... Using alternative formula for Sinch...')
%     %Sinch = eye(n2) + B^2/6 + B^4/120 + B^6/5040 + B^8/362880
% else
%     Sinch  = V*diag(sinh(diag(D))./diag(D))/V;
%     norm(inv(V)-V')
%     pause
% end 

% Sinch must be real, so we forget about imaginary parts:
Sinch = real(Sinch);

end