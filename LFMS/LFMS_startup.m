% Startup.
% Created:     12.11.2021
% Last change: 12.11.2021

close all; clear; clc;
set( 0, 'defaultAxesTickLabelInterpreter', 'latex' );
set( 0, 'defaultLegendInterpreter',        'latex' );
set( 0, 'defaultTextInterpreter',          'latex' );
set( 0, 'defaultAxesFontSize', 14 );

addpath(genpath('../LFMS'));
addpath(genpath('../altmany-export_fig'));
addpath(genpath('../Shared_Utilities'));

%--------------------------------------------------------------------------
% Define some more modern colors:
yellow = [0.894, 0.671, 0.094];
green = [0.667, 0.706, 0.118];
darkgray = [ 0.2549, 0.2549, 0.2549 ];
blue = [0.239, 0.376, 0.655];
red = [0.827, 0.341, 0.306];
orange = [0.870, 0.443, 0.137];
gray1 = [ 0.808, 0.808, 0.808 ];
gray2 = [ 0.553, 0.557, 0.557 ];
gray3 = [ 0.255, 0.255, 0.255 ];