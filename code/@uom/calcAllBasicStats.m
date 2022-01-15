function uom = calcAllBasicStats(uom)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ll = 1:numel(uom.uo)                                                   % Iterate over all the data objects
    fprintf('Getting basic stats for %d of %d:\n', ll, numel(uom.uo))
    uom.uo{ll}.calcBasicStats();                                           % Get motion errors for llth dataset, passing through parameters
    fprintf('done.\n')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%