function P = perfect_shuffle( n, m )

% function P = perfect_shuffle( n, m )
% Purpose: Returns the shuffle matrix to switch kronecker poducts
%          vec(X)  = P * vec(X'), X = n-by-m matrix 
%          and P = perfect_shuffle(m,n)
%          vec(X') = P' * vec(X)
%          so P(m,n) == P(m,n)'
% Created:     01.06.2016
% Last change: 12.07.2017


% Formulation as in book of Graham
%P = sparse(n*m,n*m);
%for i=1:n
%    for j=1:m
%        E = sparse(n,m);
%        E(i,j) = 1.;
%        P = P + kron(E,E');
%    end
%end
 
% Formulation as in Van Loan, JCAM
r = n*m;
P = [];
Ir = speye(r,r);
for i=1:m
  Ii = Ir(i:m:r,:);
  P = [P; Ii];
end
P = P';

end