function cout = pmColorMap(varargin)

% cout = pmColorMap([cSel])
%
% Funtion to create a colormap that nicely goes from one color (positive
% values) to another color (for negative values). cSel is an optional 
% input (default 'w') that specifies if the 'zero' level should be black
% ('k') or white ('w').
%
% 2021 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin > 0;    cSel = varargin{1};                                      % If provided, pass through cSel
else;             cSel = 'w';                                              % Default background color section is 'white'
end

switch lower(cSel)                                                         % Check cSel
    case 'k'                                                               % If black, leave other colors as 0
        cout = [vec(max(linspace(-1,1,255),0)), vec(max(linspace(1,-1,255),0)), 0.1*vec(linspace(0,1,255))];
    otherwise                                                              % If not black then white, 
        cout = [vec(max(linspace(-1,1,255),0)), vec(max(linspace(1,-1,255),0)), vec(abs(1*(linspace(-1,1,255))))];
        cout = 1-cout;                                                     % In this case 1-colors gives a white bg.
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% OLD CODE GRAVEYARD
% cout = [vec(max(linspace(-1,1,255),0)), vec(max(linspace(1,-1,255),0)), 0.1*ones(size(vec(linspace(-1,1,255))))];
