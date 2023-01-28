%==========================================================================
% Generates data for plotting Figure 6 in the leapfrog paper (Figure 3.6 in
% my PhD thesis). To generate the figure, run 'Driver_Figure_6_Boxplot.m'.
% Created:     20.10.2020
% Last change: 12.11.2021

%   Nov 12, 2021:
%       Added startup file.
%   Oct 20, 2020:
%       Created.
%==========================================================================

% Sets default graphics interpreter, paths and colors.
LFMS_startup;

% Global variables to store exact length and exact solution:
global Exact_Length;
global Exact_Solution;

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Set dimensions of St(n,p)
n = 12;
p = 3;
m_vector = [4:1:10,15,20:10:100];
tot_num_tests = 100;

% Fix stream of random numbers for reproducibility
s = RandStream( 'mt19937ar', 'Seed', 1 );

% % Create random Stiefel matrix Y0
X = eye( n, p );
%X = orth( rand( s, n, p ) );

% Create a random tangent vector Delta in T_{Y0}St(n,p)
distXY = 0.96*pi;
Delta_exact = distXY * GetDelta( n, p, X, s );

Exact_Length = distXY;

% Map the tg vector onto the manifold
[ Y ] = Stiefel_Exp( X, Delta_exact );

%--------------------------------------------------------------------------
% Parameters Leap-Frog and Shooting
param.tolLF = 1e-15;
param.maxiterLF = 50;
param.tolSS = 0.5e-13;
param.maxiterSS = 20;

param.verbose = 0;
%--------------------------------------------------------------------------

Q_rate_cell = cell(length(m_vector),1);

for iter_m = 1:length(m_vector)
    
    m = m_vector(iter_m);
    
    fprintf('\n----------------------------------\n')
    formatSpec = '              m = %1.0d\n';
    fprintf( formatSpec, m );
    fprintf('Random start no.: ');
    
    % Define exact solution:
    Xstar = zeros( (m-2)*n*p, 1 );
    
    for kkk=1:m-2
        
        Xstar((kkk-1)*n*p+1:kkk*n*p) = reshape( Stiefel_Exp( X, (kkk/(m-1))*Delta_exact ), [n*p,1] );
        
    end
    
    Exact_Solution = Xstar;
    
    % Legend = cell(length(m_vector),1);
    
    for test_iii = 1:tot_num_tests
        
        %     if param.verbose==1
        formatSpec = '%1.0d, ';
        fprintf( formatSpec, test_iii );
        %     end
        
        % Define initial guesses:
        Sigma_LF_0 = GetRandomStartingGuessSigmaLeapfrog( n, p, m, X, Y, param );
        
        % Build X_k:
        for k=2:m-1
            X_k((k-2)*n*p+1:(k-1)*n*p) = reshape( GetSigmakj( Sigma_LF_0, k, 1, n, p, m ), [n*p,1] );
        end
        
        % Added 16.10.2020:
        err_X_k_0 = norm( Exact_Solution - X_k', 2 );
        
        % 01.07.2020: Compute length at 0th iteration
        [ norm_error_L_omega_0, ~ ] = GetDeltaFromSigma( n, p, m, Sigma_LF_0 );
        
        % Do Leap-Frog on St(2p,p)
        [ iterLF, L_omega, norm_F_LF, norm_error_L_omega, err_X_k, Sigma_LF_k, param ] = LeapFrogStiefel_SVD( n, p, m, Sigma_LF_0, param );
        
        norm_error_L_omega = [ norm_error_L_omega_0, norm_error_L_omega ];
        err_X_k = [ err_X_k_0, err_X_k ];
        
        % This is the Q-convergence definition:
        % Q_rate(test_iii) = err_X_k(end)/err_X_k(end-1)

        % MS, 21.10.2020: Compute all the Q-rates for each k:
        Q_rate_vec_k = err_X_k(2:end)./err_X_k(1:end-1);
        
        % Save in a cell array:
        % Each Q_rate_cell{iter_m} corresponds to a different m, each
        % column of Q_rate_cell{iter_m} containes Q_rate_vec_k for a
        % different test, and each row of Q_rate_cell{iter_m} corresponds
        % to the k-th iteration
        Q_rate_cell{iter_m}(:,test_iii) = Q_rate_vec_k';
       
    end
    
end
%%
%----------------------------------------------------------------------
% SAVE DATA TO MAT-FILE
%----------------------------------------------------------------------
fileName = ['Results/Boxplot_LF_n', num2str(n), '_p', num2str(p), ...
    '_m', num2str(m_vector(1)), '_', num2str(m_vector(end)), ...
    '_tot_num_tests_', num2str(tot_num_tests), '.mat' ];
save( fileName, 'Q_rate_cell' )
fprintf('\n+--------------------------------------------------------------+\n');
fprintf('|                           Save data                          |\n');
fprintf('+--------------------------------------------------------------+\n');
fprintf('Saved data to file %s.\n', fileName );
