

% Function: phase_coherence_ts
% ----------------------------
% Calculates the phase coherence between pairs of areas in a time series matrix.
%
% Parameters:
%   - ts_phases: A matrix representing the phase values of a time series
%                 with dimensions [numAreas, numTps].
%
% Returns:
%   - iFC: A 3D matrix representing the phase coherence between pairs of areas
%           at each time point. It has dimensions [numAreas, numAreas, numTps].
%   - iFC_tril: A matrix containing the lower triangular part of iFC for each time point.
%               It has dimensions [numPairs, numTps], where numPairs is the number of
%               pairs of areas (numAreas * (numAreas - 1) / 2).
%
% Usage:
%   [iFC, iFC_tril] = phase_coherence_ts(ts_phases)
%
% Authors:
%   - Jakub Vohryzek (jakub.vohryzek@upf.edu)
%   - Yonatan Sanz-Perl (yonatan.sanz@upf.edu)
%
% Date: May 30, 2023
%
function [iFC, iFC_tril] = phase_coherence_ts(ts_phases)

    % Retrieve the dimensions of the input time series matrix
    [numAreas, numTps] = size(ts_phases);
    
    % Create an index for the lower triangular part of a matrix
    Isubdiag = find(tril(ones(numAreas), -1));
    
    % Initialising
    iFC = zeros(numAreas,numAreas,numTps);
    iFC_tril = zeros(size(Isubdiag,1),numTps);

    for t = 1:numTps
    
        iFC(:,:,t) = cos(ts_phases(:,t) - ts_phases(:,t)');
        tmp = iFC(:,:,t);
        iFC_tril(:,t) = tmp(Isubdiag);
    end
end

%% OLD VERSION
% function [iFC,iFC_tril] = phase_coherence_ts(ts_phases)
%     
%     [numAreas, numTps] = size(ts_phases);
%     Isubdiag = find(tril(ones(numAreas),-1));
% 
%     iFC = zeros(numAreas, numAreas, numTps);
%     iFC_tril = zeros(size(Isubdiag,1),numTps)
%         
%     for t = 1:numTps
%     
%         iFC(:,:,t) = cos(ts_phases(:,t) - ts_phases(:,t)');
%         tmp = iFC(:,:,t);
%         iFC_tril(:,t) = tmp(Isubdiag);
%     end
% end