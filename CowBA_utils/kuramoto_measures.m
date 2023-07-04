% Function: kuramoto_measures
% ---------------------------
% Calculates the Kuramoto order parameter and synchronization measure
% from a phase time series matrix.
%
% Parameters:
%   - phase_ts: A matrix representing the phase values of a time series
%                with dimensions [numAreas, numTps].
%
% Returns:
%   - sync: The synchronization measure, which is the mean of the Kuramoto
%           order parameter over all time points.
%   - meta: The meta-synchronization measure, which is the standard deviation
%           of the Kuramoto order parameter over all time points.
%
% Usage:
%   [sync, meta] = kuramoto_measures(phase_ts)
%
% Authors:
%   - Jakub Vohryzek (jakub.vohryzek@upf.edu)
%   - Yonatan Sanz-Perl (yonatan.sanz@upf.edu)
%
% Date: May 30, 2023
% Date: Jul 3, 2023
%
function [sync, meta,GC_proxy] = kuramoto_measures(phase_ts,amplitude_ts)

    % Retrieve the dimensions of the input time series matrix
    [numAreas, ~] = size(phase_ts);

    % Calculate the Kuramoto Order Parameter
    OP = abs(sum(exp(1i * phase_ts)) / numAreas);

    % Calculate the synchronization measure (mean of the order parameter)
    sync = mean(OP, 'omitnan');

    % Calculate the meta-synchronization measure (standard deviation of the order parameter)
    meta = std(OP, 'omitnan');
    
    % Pseudo Causality
    
    amplitude_ts = amplitude_ts/max(max(amplitude_ts));
    for seed=1:numAreas
        GC_proxy(seed) = corr2(OP(2:end),squeeze(amplitude_ts(seed,1:end-1)));
    end
    
    
end


