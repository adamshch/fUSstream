function uom = ensureAllMasks(uom, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;  
p.addParameter('redraw'      , false  );                                   % Choose to redraw the mask
parse(p,varargin{:});
p = p.Results;


for ll = 1:numel(uom.uo)                                                   % Iterate over all the data objects
    fprintf('Ensuring mask for dataset %d of %d exists...', ll, numel(uom.uo))
    uom.uo{ll}.ensureMask('redraw', p.redraw);                             % Ensure that the masks is set for the llth block
    fprintf('done.\n')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%