% Function: phases_ts
% --------------------
% Calculates the phase values of a time series matrix.
%
% Parameters:
%   - ts: A matrix representing time series data with dimensions [numAreas, numTps].
%
% Returns:
%   - ts_phases: A matrix containing the phase values of the input time series.
%                It has the same dimensions as the input [numAreas, numTps].
%
% Usage:
%   ts_phases = phases_ts(ts)
%
% Authors:
%   - Jakub Vohryzek (jakub.vohryzek@upf.edu)
%   - Yonatan Sanz-Perl (yonatan.sanz@upf.edu)
%
% Date: May 30, 2023
% Date: Jul 3, 2023
%
function [ts_phases ts_amplitude] = phases_ts(ts)

    % Retrieve the dimensions of the input time series matrix
    [numAreas, numTps] = size(ts);

    % Create a matrix of zeros with the same dimensions as the input
    ts_phases = zeros(numAreas, numTps);

    % Calculate the phase values for each area in the time series matrix
    for seed = 1:numAreas
        ts_phases(seed,:) = angle(hilbert(ts(seed,:)));
        ts_amplitude(seed,:) = abs(hilbert(ts(seed,:)))-mean(abs(hilbert(ts(seed,:))));
    end     

end


