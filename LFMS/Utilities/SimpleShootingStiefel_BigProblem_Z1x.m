function [ iter, norm_Fx_k, norm_update, Delta_k, param ] = SimpleShootingStiefel_BigProblem_Z1x( Y0, Y1, Delta_0, param )

% function [ iter, norm_Fx_k, norm_update, Delta_k, param ] = SimpleShootingStiefel_BigProblem_Z1x( Y0, Y1, Delta_0, param )
% Encapsulates all that is needed to do the simple shooting with the
% formulation using the parametrization of tangent vectors.
% Created:     22.09.2016
% Last change: 19.01.2017

% if param.verbose==1
%     disp('-----------------------------------')
%     disp('         SINGLE SHOOTING           ')
%     disp('-----------------------------------')
% end

[ n, p ] = size( Y0 );

if p==1
    % 18.01.2017: we separately handle the case of the hypersphere, i.e., St(n,1)
    Y0perp = null(Y0');
    
    Q = [ Y0, Y0perp ];
    
    d_0 = Y0perp'*Delta_0;
    
    % Initializations:
    iter = 0;
    param.flag  = true;   % we initialize the flag to true
    norm_update = param.tolSS + 1;
    norm_Fx_k   = param.tolSS + 1;
    d_k = d_0;
    
    while and( norm_update > param.tolSS, iter < param.maxiterSS )
        
        iter = iter + 1;
        
        [ Z1_d_k, J_Z1_d_k ] = GetZ1dJZ1d( Q, d_k );
        
        Fx_k = Z1_d_k - Y1(:);
        
        if cond(J_Z1_d_k) > 1e9
            if param.verbose==1
                disp( 'SINGLE SHOOTING: HUGE CONDITION NUMBER!!!' )
            end
            param.flag = false;
            break;
        end
        
        % Solve for update:
        dd_k = -J_Z1_d_k\Fx_k;
        
        d_k = d_k + dd_k;
        
        norm_update(iter) = norm( dd_k );
        norm_Fx_k(iter)   = norm( Fx_k );
        
    end
    
    if iter==param.maxiterSS
        if param.verbose==1
            disp( 'SINGLE SHOOTING: REACHED MAXIMUM NUMBER OF ITERATIONS!!!' )
        end
        param.flag = false;
        Delta_k = 0;
        return;
    end
    
    Delta_k = Y0perp*d_k;
    
else
    
    % The bases for the space Skew(p)
    E = GetBasisSkewp( p );
    
    % Dimension of Skew(p)
    dimSkewp = (1/2)*p*(p-1);
    
    Y0perp = null(Y0');
    
    Q = [ Y0, Y0perp ];
    
    % 16.01.2017: I removed the initialization of Delta_0 here and set it as an
    %             input of this function.
    
    % "Condense" into the vector of p(n-(1/2)(p+1)) parameters x_0
    Omega_0 = Y0'*Delta_0;
    s_0 = E\Omega_0(:);
    
    D_0 = Y0perp'*Delta_0;
    d_0 = D_0(:);
    
    x_0 = [ s_0; d_0 ];
    
    % Initializations:
    iter = 0;
    param.flag  = true;   % we initialize the flag to true
    norm_update = param.tolSS + 1;
    norm_Fx_k   = param.tolSS + 1;
    %rate        = param.tolSS + 1;
    x_k = x_0;
    
    while and( norm_update > param.tolSS, iter < param.maxiterSS )
        
        iter = iter + 1;
        
        % The Jacobian
        J_Z1_x = GetJZ1xAnalytic( E, Q, x_k );
        
        % J_Z1_x_plus = pinv(J_Z1_x);   % only for analysis purposes
        
        % The function F(Delta)
        Fx_k = GetZ1x( E, Q, x_k ) - Y1(:);
        
        % norm( reshape(Fx_k,[n,p]))
        % norm( Fx_k )
        % normJ_Z1_x_plus=norm(J_Z1_x_plus)
        % pause
        % min(svd(J_Z1_x_plus))
        % norm_Prod = norm( J_Z1_x_plus * Fx_k)
        
        % 18.12.2016: We put a check on the condition number of the Jacobian
        if cond(J_Z1_x) > 1e9
            if param.verbose==1
                disp( 'SINGLE SHOOTING: HUGE CONDITION NUMBER!!!' )
            end
            param.flag = false;
            break;
        end
        
        % Solve for update:
        dx_k = -J_Z1_x\Fx_k;
        
        %check_pinv = - J_Z1_x_plus*Fx_k - dx_k
        
        x_k = x_k + dx_k;
        
        norm_update(iter) = norm( dx_k );
        norm_Fx_k(iter)   = norm( Fx_k );

        % 23.01.2017: stop if norm_Fx_k is increasing
%         if iter==10
%             if norm_Fx_k(iter) > norm_Fx_k(1)
%                 if param.verbose==1
%                     disp( 'SINGLE SHOOTING: norm_Fx_k is increasing!!!' )
%                 end
%                 param.flag = false;
%                 break;
%             end
%         end
        
        %     if iter~=1
        %         rate(iter) = norm_update(iter)/norm_update(iter-1);
        %     end
        
        %     % 19.09.2016 -----------------------------
        %     fprintf( 'Norm of A:       %1.6f \n', norm( A ) )
        %     fprintf( 'Norm of matexpA: %1.6f \n', norm( matexpA ) )
        %     fprintf( 'Norm of KA:      %1.6f \n', norm( KA ) )
        %     fprintf( 'Norm of J:       %1.6f \n', norm( J ) )
        %     %-----------------------------------------
        %     if param.verbose==1
        %         formatSpec = 'iter %3.2d   cond(J_Z1_x) = %5.3f   norm_update = %10.8f   norm_Fx_k = %10.8f   rate = %10.8f \n';
        %         fprintf( formatSpec, iter, cond(J_Z1_x), norm_update(iter), norm_Fx_k(iter), rate(iter) )
        %     end
        
    end
    
    if iter==param.maxiterSS
        % 18.12.2016: if we reach the maximum number of iterations allowed for
        % Single Shooting, we set the flag to 'false'
        if param.verbose==1
            disp( 'SINGLE SHOOTING: REACHED MAXIMUM NUMBER OF ITERATIONS!!!' )
        end
        param.flag = false;
        Delta_k = 0;
        return;
    end
    
    % Build Delta_k from the parameters s and d:
    % Extraction of s and d from x_k:
    s = x_k(1:dimSkewp);
    d = x_k(dimSkewp+1:end);
    
    vecOmega = E*s;
    
    % Matricization of vecOmega
    Omega = reshape( vecOmega, [p,p] );
    
    % Matricization of d
    D = reshape( d, [n-p,p] );
    
    Delta_k = Y0*Omega + Y0perp*D;
end
end


function [ Z1_d, J_Z1_d ] = GetZ1dJZ1d( Q, d )

n = size(Q,1);

if n==2
    % We are in the very special case St(2,1):
    Z1_d   = Q * [  cos(d); sin(d) ];
    J_Z1_d = Q * [ -sin(d); cos(d) ];
else
    %--------------------------------------------------------------------------
    % Initialize matrices that depend only on n and p
    IpZeros = eye( n, 1 );
    % Transformation matrix from blockwise vectorization to ordinary vectorization
    H = GetH( n, 1 );
    %--------------------------------------------------------------------------
    
    % Build matrix h(d)
    hd = [ 0,            -d';
        d,    zeros( n-1 ) ];
    
    % matexphx
    matexphd_k = expm( hd );
    
    % Geodesics
    Z1_d = kron( [1, zeros(1,n-1)], Q ) * matexphd_k(:);
    
    % Get Kronecker representation of the Fréchet derivative of the matrix
    % exponential
    Khd = GetKA( hd );
    
    % Build Jacobian of h(d) wrt d
    J_hd = [   zeros(1,n-1);
        eye(n-1);
        -eye(n-1);
        zeros((n-1)^2, n-1) ];
    
    % The Jacobian of Z1 wrt d
    J_Z1_d = kron( IpZeros', Q ) * Khd * H * J_hd;
end

end