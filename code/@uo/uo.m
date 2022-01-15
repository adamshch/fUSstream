classdef uo < handle
    
%    This file defines the functional ultrasound object uo
%
% 2020 - Adam Charles 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties                                                             % Set up UO object properties
        
        movie                                                              % 3D movie object or empty if meta has pointer to matfile (keeps movie on disk)
        eventVec                                                           % Vector of events
        frameTimes                                                         % Times at which frames occured
        frameSize                                                          % Size of a single frame (2x1 vector: depth x AP length)
        scanArea                                                           % Area that is scanned (depth x AP length)
        scanRes                                                            % Scanning resolution (depth x AP length)
        movieLen                                                           % Duration of movie in # of frames
        dt                                                                 % Movie Sampling rate
        mask                                                               % Mask of where the brain is: must be drawn or provided. 
        stats  = struct('medImg',[],'meanImg',[],'maxImg',[],...
                                                 'minImg',[],'varImg',[]); % Place to store basic statistics about the fu movie
        meta   = struct('fileName',[],'dataPath',[],'savePath',[],...
                        'subjectID',[],'movVarName',[],'loadMov',[],...
                                                         'movMatFile',[]); % Meta substructure keeps extra info related to the experiment etc. 
        errs   = struct('burst', [], 'motion', []);                        % Error statistics
        post   = struct('movDenoised', [],'movMotionCorrected', []);       % post-processing results
        PCAout = struct();                                                 % PCA output
        GraFTout = struct();                                               % GraFT output
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Below are the set of methods associated with this class

    methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Main function
        function uo = uo(movieSpec, varargin)
            p = inputParser;                                               % Set up an object to parse all of the various inputs
            p.addParameter('nFrames'      , []);                           % Can opt not to load the entire movie
            p.addParameter('loadToRAM'    , false);                        % Choose if to load the entire data to RAM or to just link to the MAT file
            p.addParameter('scanArea'     , []);                           % Need to input either the scan area...
            p.addParameter('scanRes'      , []);                           % ...or scan resolution
            p.addParameter('dt'           , []);                           % Option to give the temporal scan resolution...
            p.addParameter('frameTimes'   , []);                           % ...and/or the times of frames.
            p.addParameter('preCalcStats' , false);                        % ...and/or the times of frames.
            p.addParameter('events'       , []);                           % Option to input an event vector. If specified as a non-empty array will override the loaded event vector
            p.addParameter('drawMask'     , true);                         % Option to draw a mask over the brain
            parse(p,varargin{:});
            p = p.Results;
            
            uo.meta.loadMov = p.loadToRAM;                                 % Save in the object whether the movie should ever be loaded. 
            
            if isempty(movieSpec)                                          % If nothing is given...
            elseif isnumeric(movieSpec)                                    % If a 3D array (the actual movie) was passed in...
                assert(ndims(movieSpec)==3)                                % ...Ensure that the movie is a 3D array
                movSz         = size(movieSpec);                           % ...Get the size of the movie
                uo.movie      = movieSpec;                                 % ...Include movie in the object
                uo.frameSize  = movSz(1:2);                                % ...Get the size of each frame
                uo.movieLen   = movSz(3);                                  % ...Get the movie duration
                uo.eventVec   = p.events;                                  % ...Pass through the event vector
                uo.dt         = p.dt;                                      % ...Pass through the temporal sampling rate
                uo.frameTimes = p.frameTimes;                              % ...Pass through the frame times
                uo            = checkFrameTimes(uo, p);                    % Check the parameters for frame rate and time
                uo            = makeDataAndSavePath(uo, []);               % ...Create the different directories needed for saving
                uo.meta.movMatFile = [];
                uo.meta.loadMov    = true;
            elseif isstring(movieSpec)||ischar(movieSpec)                  % If a string (path of the movie) was passed in...
                assert(isfile(movieSpec));                                 % ...Make sure the file actually exists
                assert(strcmp(movieSpec(end-3:end),'.mat'));               % ...For now only mat files are permitted
                uo = loadUltrasoundDataFromMat(uo, movieSpec, p);          % ...Load the movie and events vector
                if p.loadToRAM
                    movSz        = size(uo.movie);                         % ...Get the size of the movie
                    uo.frameSize = movSz(1:2);                             % ...Get the size of each frame
                    uo.movieLen  = movSz(3);                               % ...Get the movie duration
                    uo.eventVec  = checkEventVec(uo.eventVec, p.events);
                    uo           = checkFrameTimes(uo, p);
                    uo           = makeDataAndSavePath(uo, movieSpec);     % ...Extract directory and file information, and create
                else
                    movSz        = size(...
                                  uo.meta.movMatFile.(uo.meta.movVarName));% ...Get the size of the movie
                    uo.frameSize = movSz(1:2);                             % ...Get the size of each frame
                    uo.movieLen  = movSz(3);                               % ...Get the movie duration
                    uo.eventVec  = checkEventVec(uo.eventVec, p.events);
                    uo           = checkFrameTimes(uo, p);
                    uo           = makeDataAndSavePath(uo, movieSpec);     % ...Extract directory and file information, and create
                end
            else
                error('Cannot recognize movieSpec: invalid entry!')
            end
            
            if p.preCalcStats
                uo = uo.calcBasicStats();
            end
            uo = setResolution(uo, p, uo.frameSize);
            
        end
        function uo = makeDataAndSavePath(uo, pathFull)
            if isempty(pathFull)
                uo.meta.dataPath  = [];                                    % The current data path is unknown! 
                uo.meta.fileName  = [];                                    % The current data file name is unknown! 
                saveTmp           = pwd;                                   % Create a potential save name that follows the basic convention
                saveTmp           = [saveTmp,'/',datestr(now,'yyyymmdd'),...
                                                             'fuResults']; % Make a subfolder based on the date
                uo.meta.savePath  = saveTmp;                               % save the subfolder to the object
                if ~isfolder(saveTmp)
                    mkdir(saveTmp)
                end
            else
                [uo.meta.dataPath, uo.meta.fileName] = fileparts(pathFull);% Extract the data path and file name
                saveTmp = strrep(uo.meta.dataPath,'Data','Analysis');      % Create a potential save name that follows the basic convention
                if isfolder(saveTmp)
                    uo.meta.savePath = saveTmp;                            % Set the directory which all the results are saved to
                end
            end
        end
        function uo = setResolution(uo, p, frameSize)
            if isempty(p.scanArea)&&isempty(p.scanRes)
                warning('Neither the scan area or scan resultion are set! Defaulting to a typical 15x16mm scan area')
                uo.scanArea = [15,16];
                uo.scanRes  = [15,16]./frameSize(:)';
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      This function calculates a set of basic statistics from the movie.
       varargout = getBasicStats(uo, varargin)
       uo = calcBasicStats(uo, varargin)
       getTrace = makeGetTraceFunction(uo,p)
       getBlock = makeGetBlockFunction(uo,p)
       getMovie = makeGetMovieFunction(uo,p)
       uo = drawROI(uo, varargin)
       uo = ensureMask(uo, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      This function finds burst errors (frames with excessively high
%      error rates):
       varargout = findBurstFrames(uo, varargin);                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      This function finds motion correction errors (frames where the
%      motion correction failed):
       varargout = findMotionCorrectionError(uo, varargin)
       uo = correctResidualMotion(uo, varargin)
       motCorr = correlateMotionWithData(uo, varargin)
       motTest = motionHypothesisTest(uo, varargin)
       motTest = compareMotionHypothesisTest(uo, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      This function plots examples of found burst errors. Must run
%      'findBurstFrames' above first (or provide a vector of logical
%      values for uo.
       plotBurstErrors(uo, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       displayMovie(uo, varargin)                                          % Display the movie using MovieSlider
       displayExampleFrame(uo, varargin)                                   % Plot an example frame
       varargout = displayExampleTraces(uo, varargin)                      % Plot example time traces
       displayMotionError(uo, varargin)
       writeToAVI(uo, video_name, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       denoiseTraces(uo, varargin)                                         % Wavelet denoising of traces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       fuPCA(uo, varargin)                                             % Run PCA on the data
       % uo = isomapProj(uo, varargin)
       % uo = graphMap(uo, varargin)
       % uo = correlationMap(uo, varargin)
       % uo = clusterTraces(uo, varargin)
       % EventShape        = event2timeseries(uo, varargin)
       % [eventframe, idx] = event2frame(uo, events, varargin)
       eSTH = eventSTH(uo, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       varargout = computeCorrsWithBaseline(uo, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % TFout = doesDenoiseExist(uo)
       % uo = ensureDenoisedData(uo)
       % TFout = doesMotionCorrectedExist(uo)
       % TFout = isMatFileBased(uo)
       % uo = ensureMedian(uo)
       % uo = ensureMotionErrsComputed(uo)
       % motCorrSize = getMotionCorrectedSize(uo)
       % [motDep, motLat, maxMot] = getMotionVectors(uo)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uo = loadUltrasoundDataFromMat(uo, fileName, p, varargin)

% mov = loadUltrasoundDataFromMat(fileName, varargin)
%
% Function to load up an ultrasound dataset
%
% 2019 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin > 3
    nameToLoad = varargin{1};
else
    nameToLoad = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data


m                  = matfile(fileName);                                    % Create mat file object to load data from
uo.meta.movMatFile = m;                                                    % Save the mat file object that the data is loaded from

if isempty(nameToLoad)
    fNames    = fieldnames(m);                                             % Get all the field names
    fieldSize = cell(numel(fNames),2);                                     % Initialize a cell array to store sizes of all mat variables
    for ll = 1:numel(fNames)
        fieldSize{ll,1} = size(m.(fNames{ll}));                            % Get the size of all objects in the mat file
        fieldSize{ll,2} = numel(fieldSize{ll,1});                          % Get the dimensions of all objects in the mat file
    end
    idx3D      = cell2mat(fieldSize(:,2)) == 3;                            % Find which arrays have 3 dimensions
    assert(any(idx3D))                                                     % Make sure there is at least one 3D array
    allSzs     = cellfun(@prod, fieldSize(:,1));                           % get the number of elements for each variable
    max3D      = (allSzs == max(allSzs(idx3D)))&idx3D;                     % Pick out the 3D array with the most elements (most likely to be the video)
    uo.meta.movVarName = fNames{max3D};                                    % Pull out the name of the biggest 3D object
    nameToLoad = uo.meta.movVarName;
end

if isempty(p.nFrames); p.nFrames = fieldSize{max3D,1}(3);  end             % If nFrames is empty then load all the frames

if p.loadToRAM
    if numel(p.nFrames) == 1
        uo.movie = m.(nameToLoad)(:,:,1:p.nFrames);
    elseif numel(p.nFrames) == 2
        uo.movie = m.(nameToLoad)(:,:,min(p.nFrames):max(p.nFrames));
    elseif numel(p.nFrames) > 2
        uo.movie = m.(nameToLoad)(:,:,p.nFrames);
    else
        error('Invalid frame selection!')
    end
else
    uo.movie = [];
end

if any(strcmp(fNames,'events'));    uo.eventVec = m.events;
else;                               uo.eventVec = [];
end

if any(strcmp(fNames,'frames'));    uo.frameTimes = m.frames;
else;                               uo.frameTimes = [];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function eventVec  = checkEventVec(eventVec, pEvents)

if isstring(pEvents)
    assert(isfile(pEvents))                                                % Make sure the file actually exists
    assert(strcmp(pEvents(end-3:end),'.mat'))                              % ...For now only mat files are permitted
    m       = matfile(pEvents);                                            % Access the mat file
    pEvents = m.events;                                                    % Get the numeric array of the provided
end

if ~isempty(pEvents); eventVec = pEvents; end                              % Pass through the event vector

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function uo  = checkFrameTimes(uo, p)

if isempty(uo.frameTimes)&&isempty(uo.dt)                                  % If nothing is provided...
    if (~isempty(p.dt))||(~isempty(p.frameTimes))
        uo.dt         = p.dt;                                              % Set pre-defined sampling rate
        uo.frameTimes = p.frameTimes;                                      % Set pre-defined frame times
        uo            = checkFrameTimes(uo, p);                            % Re-check parameters
    else
        warning('No timing information given. Assuming a 2Hz frame-rate')
        uo.dt         = 0.5;                                               % Set default sampling rate
        uo            = checkFrameTimes(uo, p);                            % Re-check parameters
    end
elseif (~isempty(uo.frameTimes))&&(isempty(uo.dt))                         % If only frameTimes is provided...
    uo.dt = mean(diff(uo.frameTimes));
    if var(diff(uo.frameTimes))>uo.dt
        warning('Variance of time intervals is greater than the average time interval! Take dt estimate with caution!')
    end
elseif (isempty(uo.frameTimes))&&(~isempty(uo.dt))                         % If only dt is provided...
    uo.frameTimes = (0:(uo.movieLen-1))/uo.dt;                             % ...Create time vector from temporal sampling interval
elseif (~isempty(uo.frameTimes))&&(~isempty(uo.dt))                        % If both are provided...
    assert(numel((0:(uo.movieLen-1))/uo.dt)==numel(uo.frameTimes))         % ... Make sure the sizes fit
else;  error('Something is very wrong if you get this error...')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
