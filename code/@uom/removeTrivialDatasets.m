function uom = removeTrivialDatasets(uom)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

keepDatasets = true(numel(uom.uo),1);

for ll = 1:numel(uom.uo)                                                   % Iterate over all the data objects
    if norm(uom.uo{ll}.stats.medImg(:))<eps
        keepDatasets(ll) = false;
    end
end

uom.uo = uom.uo(keepDatasets);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%