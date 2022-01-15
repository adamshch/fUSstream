function uo = ensureMask(uo, varargin)

% function uo = ensureMask(uo)
%
% Function to ensure that a mask is specified
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;  
p.addParameter('redraw'      , false  );                                   % Choose to redraw the mask
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(uo.mask)
    warning('Mask not defined. Please specify mask...\n')
    uo.drawROI('redraw', p.redraw);
else
    if p.redraw
        fprintf('Mask exists but redrawing anyway...\n')
        uo.drawROI('redraw', p.redraw);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%