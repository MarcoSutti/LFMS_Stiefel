function [ Delta_0 ] = GetStartingGuessDelta( Y0, Y1 )

% function [ Delta_0 ] = GetStartingGuessDelta( Y0, Y1 )
% Purpose: Returns the initial guess Delta_0 for Single Shooting.
% Created:     16.01.2017
% Last change: 16.01.2017

% Initial guess for Delta:  (26.08.2016)
DeltaFO = Y1 - Y0;
alpha = norm( DeltaFO );
ProjDeltaFO = ProjTgSpaceStiefel( Y0, DeltaFO );
beta = norm( ProjDeltaFO );
Delta_0 = (alpha/beta) * ProjDeltaFO;

end