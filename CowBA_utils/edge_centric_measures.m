% Function: edge centric measures
% ---------------------------
% Calculates the Edge centric mestastability
% from a  time series matrix.
%
% Parameters:
%   - ts filtered: A matrix representing the  values of a time series
%                with dimensions [numAreas, numTps].
%
% Returns:
%   - EdgeMeta: The synchronization measure, which is the mean of the Kuramoto
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
function [EdgeMeta] = edge_centric_measures(ts)

    % Retrieve the dimensions of the input time series matrix
    [numAreas, Tmax] = size(ts);

    for t=1:Tmax
        EdFC=ts(:,t)*ts(:,t)';
        EdgesEmp(:,t)=EdFC(:);
    end
    EdgeMeta = std(EdgesEmp(:),'omitnan');
  
    
end


