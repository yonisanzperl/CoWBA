% Function: filter_ts
% -------------------
% Applies a filter to a time series matrix, handling NaN values in the input.
%
% Parameters:
%   - ts: A matrix representing time series data with dimensions [numAreas, numTps].
%   - bfilt: Numerator coefficients of the filter.
%   - afilt: Denominator coefficients of the filter.
%
% Returns:
%   - ts_filtered: A matrix containing the filtered time series data.
%                  It has the same dimensions as the input [numAreas, numTps].
%
% Usage:
%   ts_filtered = filter_ts(ts, bfilt, afilt)
%
% Authors:
%   - Jakub Vohryzek (jakub.vohryzek@upf.edu)
%   - Yonatan Sanz-Perl (yonatan.sanz@upf.edu)
%
% Date: May 30, 2023
%
function [ts_filtered] = filter_ts(ts, bfilt, afilt)

    % Retrieve the dimensions of the input time series matrix
    [numAreas, numTps] = size(ts);

    % Create a matrix of zeros with the same dimensions as the input
    ts_filtered = zeros(numAreas, numTps);

    % Apply the filter to each area in the time series matrix
    for seed = 1:numAreas

        % Check if there are any NaN values in the time series
        if sum(isnan(ts(seed,:))) < 1
            % If there are no NaN values, apply the filter
            ts_filtered(seed,:) = filtfilt(bfilt, afilt, ts(seed,:));
        else
            % If there are NaN values, keep the original time series without filtering
            ts_filtered(seed,:) = ts(seed,:);
        end
    end     

end


