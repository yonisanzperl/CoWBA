% Function: LZc_measure
% ----------------------
% Calculates the normalized Lempel-Ziv Complexity (LZc) measure for each area
% in a time series matrix.
%
% Parameters:
%   - ts: A matrix representing the time series with dimensions [numAreas, numTps].
%
% Returns:
%   - C_norm: A vector containing the normalized LZc measure for each area.
%
% Usage:
%   C_norm = LZc_measure(ts)
%
% Authors:
%   - Jakub Vohryzek (jakub.vohryzek@example.com)
%   - Yonatan Sanz-Perl (yonatan.sanzperl@example.com)
%
% Date: May 30, 2023
% Date: Jul 3, 2023


function [C_norm,Cmean] = LZc_measure(ts)

    % Retrieve the dimensions of the input time series matrix
    [numAreas, ~] = size(ts);

    % Preallocate C_norm variable
    C_norm = zeros(1,numAreas);
    
    % Calculate Lempel-Ziv Complexity
    for i = 1:numAreas
        tmp = ts(i, :);

        % Apply thresholding: setting all values above 0 to 1, else to 0
        tmp = tmp > 0;

        [C_norm(i), ~, ~] = calc_lz_complexity(tmp, 'exhaustive', true);
    end
    Cmean= mean(C_norm);
end
