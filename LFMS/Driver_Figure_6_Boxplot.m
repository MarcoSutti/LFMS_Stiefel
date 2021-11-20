%==========================================================================
% Generates Figure 6 in the leapfrog paper (Figure 3.6 in my PhD thesis).
% For post-processing Q_rate_matrix.
% To generate the data, run 'Driver_Leapfrog_random_starts_2020.m'.
% Created:     21.10.2020
% Last change: 12.11.2021

%   Nov 12, 2021:
%       Added startup file.
%   Jun 4, 2021:
%       Modified name of file and description.
%   Oct 21, 2020:
%       Created.
%==========================================================================

% Sets default graphics interpreter, paths and colors.
LFMS_startup;

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Set dimensions of St(n,p)
n = 12;
p = 3;
m_vector = [4:1:10,15,20:10:100];
tot_num_tests = 100;


% Decide which entries to extract from m_vector:
idx_m_vector = [1:2:7, 8:length(m_vector)]


m_vec_reduced = m_vector(idx_m_vector)

%--------------------------------------------------------------------------
length_m_vec = length(m_vec_reduced);
%--------------------------------------------------------------------------
% Open the file
% This loads the cell array Q_rate_cell{}
load('Results/Boxplot_LF_n12_p3_m4_100_tot_num_tests_100.mat')

% Each Q_rate_cell{iter_m} corresponds to a different m, each column of
% Q_rate_cell{iter_m} containes Q_rate_vec_k for a different test, and each
% row of Q_rate_cell{iter_m} corresponds to the k-th iteration.

worst_Q_rate_matrix = zeros(tot_num_tests,length_m_vec);
median_Q_rate_matrix = zeros(tot_num_tests,length_m_vec);
Q_rate_1_matrix = zeros(tot_num_tests,length_m_vec);
Q_rate_end_matrix = zeros(tot_num_tests,length_m_vec);


for iter_m = 1:length_m_vec

    for test_iii = 1:tot_num_tests
        Q_rate_vec_k = Q_rate_cell{iter_m}(:,test_iii);
        
        % Compute the worst:
        worst_Q_rate_k(test_iii) = max( Q_rate_vec_k );
        
        % Compute the median:
        median_Q_rate_k(test_iii) = median( Q_rate_vec_k );
        
        % Compute the first:
        Q_rate_1(test_iii) = Q_rate_vec_k(1);
        
        % Compute the last:
        Q_rate_end(test_iii) = Q_rate_vec_k(end);
    
    end    
    
    % For each m, store the worst, the median, the first and the last
    % Q-rate for all the tests. Each row corresponds to a different tests,
    % each column to a different m.
    worst_Q_rate_matrix( :, iter_m ) = worst_Q_rate_k;
    median_Q_rate_matrix( :, iter_m ) = median_Q_rate_k;
    Q_rate_1_matrix( :, iter_m ) = Q_rate_1;
    Q_rate_end_matrix( :, iter_m ) = Q_rate_end;

end

%--------------------------------------------------------------------------
figure(1);

boxplot( worst_Q_rate_matrix, 'BoxStyle', 'outline', 'MedianStyle', 'line', 'Widths', 0.5 ) %, 'Color', orange )
handle_outliers = findobj(gca,'tag','Outliers');   % Get handles for outlier lines.
handle_box = findobj(gca,'tag','Box');        % Get handles for the Box
handle_median = findobj(gca,'tag','Median');    % Get handles for the median line
handle_uw = findobj(gca,'tag','Upper Whisker');    % Get handles for the upper whisker
handle_lw = findobj(gca,'tag','Lower Whisker');    % Get handles for the lower whisker
handle_uav = findobj(gca,'tag','Upper Adjacent Value');
handle_lav = findobj(gca,'tag','Lower Adjacent Value');

% Modify colors and line widths of the various component of the boxplot
set( handle_outliers, 'Marker', 'o', 'MarkerEdgeColor', gray2, 'LineWidth', 1.5 );
set( handle_box, 'Color', blue, 'LineWidth', 2 );    % Change color for box
set( handle_median, 'Color', red, 'LineWidth', 2 );
set( handle_uw, 'Color', gray2, 'LineWidth', 1.5 );
set( handle_lw, 'Color', gray2, 'LineWidth', 1.5 );
set( handle_uav, 'Color', gray2, 'LineWidth', 1.5 );
set( handle_lav, 'Color', gray2, 'LineWidth', 1.5 );
set(gca,'box','off')


hold on
plot( [ 0 100 ], [1 1], '--', 'Color', gray2, 'LineWidth', 2 )

xticklabels( m_vec_reduced )

% Legend
legend( [ handle_median(1), handle_box(1), handle_outliers(1) ], ...
    {'Median', '25\%-75\%', 'outliers'}, 'FontSize', 14, 'Location', 'SE' );


xlabel( 'number of points $m$' )
ylabel('$\mu_{k}$')

set(gca, 'TickLabelInterpreter', 'latex');

fileName = ['Plots/Boxplot_LF_n', num2str(n), '_p', num2str(p), ...
    '_m', num2str(m_vec_reduced(1)), '_', num2str(m_vec_reduced(end)), ...
    '_tot_num_tests_', num2str(tot_num_tests) ];

saveas( gcf, fileName, 'epsc' )
fprintf('-------------------------------------------------------------------------\n');
fprintf('Saved graph to file %s.eps.\n', fileName);
fprintf('-------------------------------------------------------------------------');
% export_fig Plots/Boxplot_Leapfrog_m_10_100.pdf -pdf -cmyk -transparent;

worst_of_the_best_Q_rate = max(Q_rate_1_matrix)
worst_of_the_worst_Q_rate = max(worst_Q_rate_matrix)
median_of_the_worst_Q_rate = median(worst_Q_rate_matrix)

%--------------------------------------------------------------------------
figure(2)

% Plot Q-convergence factor versus m
handle_array = plot( m_vec_reduced, worst_of_the_worst_Q_rate, 'd-', 'Color', ...
    blue, 'LineWidth', 2, 'MarkerEdgeColor', blue, 'MarkerFaceColor', ...,
    blue, 'MarkerSize', 8 );
hold on
plot( [m_vec_reduced(1) m_vec_reduced(end)], [1 1], '--', 'Color', gray2, 'LineWidth', 2 )
grid on
xlim( [m_vec_reduced(1) m_vec_reduced(end)] )
ylim( [0.865, 1.01])
xlabel( 'number of points $m$' )
ylabel('$\mu_{k}$')

drawnow;
% LineWidth of the MarkerEdge:
myMarkerLineWidth = 0.5;
handle_array.MarkerHandle.LineWidth = myMarkerLineWidth;

export_fig Plots/Convg_Qrate_vs_m.pdf -pdf -cmyk -transparent;
