function EventShape = event2timeseries(fuObj, varargin)

%
%
%
%
% 2020 - Ahmed El-Hady & Adam Charles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of th e various inputs
% p.addParameter('frames' , (0:(fuObj.movieLen-1)));                         % Can opt not to analyze all the frames and instead select which frames align with the events vector
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[eventframe, idx] = fuObj.event2frame();
EventShape        = zeros(length(fuObj.frameTimes),1);

for ne = 1:size(eventframe,1)
    EventShape(idx(ne,1):idx(ne,2)) = 1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%