function  [eventframe, idx] = event2frame(fuObj, varargin)

%
%
% Given the event and frame time stamps is seconds (or any other unfiromly 
% spaced time unit) event2frame returns the frame index at which the event 
% happened. Given that the sampling rate of fUSi imaging is around 2Hz, 
% the precision of the event is less than 250ms.
%
% 2020 - Ahmed El-Hady & Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of th e various inputs
% p.addParameter('frames'  , 1:fuObj.movieLen);                              % Can opt not to analyze all the frames and instead select which frames align with the events vector
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations

ev = fuObj.eventVec(:);                                                    % Vectorize the events matrix
le = size(fuObj.eventVec, 1);                                              % Get number of events

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main block

if min(size(fuObj.eventVec)) == 2
    idx = NaN(le,2);
    for ne = 1:le
        [~, idx(ne,1)] = min(abs(ev(2*ne-1) - fuObj.frameTimes));                  % checking the closest frame
        [~, idx(ne,2)] = min(abs(ev(2*ne) - fuObj.frameTimes));                    % checking the closest frame
    end
    eventframe = fuObj.frameTimes(idx);
    
elseif min(size(fuObj.eventVec)) == 1
    idx = NaN(le,1);
    for ne = 1:le
        [~, idx(ne)] = min(abs(ev(ne) - fuObj.frameTimes));                        % checking the closest frame
    end
    eventframe = fuObj.frameTimes(idx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%