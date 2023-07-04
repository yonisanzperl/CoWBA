% Function: demean_detrend_ts
% ---------------------------
% Removes the mean and trends from a time series matrix.
%
% Parameters:
%   - ts: A matrix representing time series data with dimensions [numAreas, numTps].
%
% Returns:
%   - ts_processed: A matrix containing the demeaned and detrended time series data.
%                   It has the same dimensions as the input [numAreas, numTps].
%
% Usage:
%   ts_processed = demean_detrend_ts(ts)
%
% Authors:
%   - Jakub Vohryzek (jakub.vohryzek@upf.edu)
%   - Yonatan Sanz-Perl (yonatan.sanz@upf.edu)
%
% Date: May 30, 2023
%
function [ts_processed] = demean_detrend_ts(ts)

    % Retrieve the dimensions of the input time series matrix
    [numAreas, numTps] = size(ts);

    % Create a matrix of zeros with the same dimensions as the input
    ts_processed = zeros(numAreas, numTps);

    % Remove the mean and trends for each area in the time series matrix
    for seed = 1:numAreas
        ts_processed(seed,:) = detrend(ts(seed,:) - mean(ts(seed,:), 'omitnan'));
    end     

end

