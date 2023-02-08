function varargout = correlateToEvent(uo, varargin)

% varargout = correlateToEvent(uo, varargin)
%
% 
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('numClusters'   , 20);                                      % Number of clusters
p.addParameter('motionCorrect' , false   );                                % Select if the motion corrected data should be used
p.addParameter('denoised'      , false   );                                % Choose to display denoised movie
p.addParameter('plotOpt'       , true    );                                % Choose to display denoised movie
p.addParameter('useMask'       , true    );
p.addParameter('window'        , [3,5]   );                                % 1 or 2 parameter window (either symmetric about the point or pre/post point
% p.addParameter('useMask'       , true    );

parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

EventShape        = uo.event2timeseries(uo, varargin);
[eventframe, idx] = uo.event2frame(uo, events, varargin);



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
