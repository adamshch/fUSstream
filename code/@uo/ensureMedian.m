function uo = ensureMedian(uo)

% function uo = ensureMedian(uo)
%
% Function to ensure that the median image is pre-computed
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isequal(size(uo.stats.medImg),uo.frameSize)
    warning('Median not pre-computed. Computing basic sats...')
    uo.calcBasicStats();
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%