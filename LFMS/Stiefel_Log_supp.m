function [Delta, k, conv_hist, norm_logV0] = Stiefel_Log_supp(U0, U1, tau)
%-------------------------------------------------------------
%@author: Ralf Zimmermann, IMADA, SDU Odense
%
% Input arguments
% U0, U1 : points on St(n,p)
%    tau : convergence threshold
% Output arguments
%  Delta : Log^{St}_U0(U1),
%          i.e., tangent vector such that Exp^St_U0(Delta) = U1
%      k : iteration count upon convergence
% supplementary output
%  conv_hist : convergence history
% norm_logV0 : norm of matrix log of first iterate V0
%-------------------------------------------------------------

% get dimensions
[n,p] = size(U0);
% store convergence history
conv_hist = [0];

% step 1
M = U0'*U1;
% step 2
[Q,N] = qr(U1 - U0*M,0);   % thin qr of normal component of U1
% step 3
[V, ~] = qr([M;N]);                   % orthogonal completion

% "Procrustes preprocessing"
[D,S,R]      = svd(V(p+1:2*p,p+1:2*p));
V(:,p+1:2*p) = V(:,p+1:2*p)*(R*D');
V            = [[M;N], V(:,p+1:2*p)];  %          |M  X0|
                                       % now, V = |N  Y0|
% just for the record
norm_logV0 = norm(logm(V),2);

% step 4: FOR-Loop
for k = 1:10000
    % step 5
    [LV, exitflag] = logm(V);
                                 % standard matrix logarithm
                                 %             |Ak -Bk'|
                                 % now, LV =   |Bk  Ck |
    C = LV(p+1:2*p, p+1:2*p);    % lower (pxp)-diagonal block

    % steps 6 - 8: convergence check
    normC = norm(C, 2);
    conv_hist(k) = normC;
    if normC<tau;
        disp(['Stiefel log converged after ', num2str(k),...
              ' iterations.']);
        break;
    end
    
    % step 9
    Phi = expm(-C);              % standard matrix exponential
    % step 10
    V(:,p+1:2*p) = V(:,p+1:2*p)*Phi;   % update last p columns
end

% prepare output                         |A  -B'|
% upon convergence, we have  logm(V) =   |B   0 | = LV
%     A = LV(1:p,1:p);     B = LV(p+1:2*p, 1:p)
% Delta = U0*A+Q*B
Delta = U0*LV(1:p,1:p) + Q*LV(p+1:2*p, 1:p);
return;
end