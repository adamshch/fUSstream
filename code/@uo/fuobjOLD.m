% classdef fuobjOLD < handle
    
    
%    This file defines the functional ultrasound object fuobj
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        
        movie      % 3D movie object
        eventVec   % Vector of events
        frameTimes % Times at which frames occured
        dataPath   % Path to the data
        savePath   % Path where all results are saved
        fileName   % Name of file where data exists
        subjectID  % Name of subject the movie was taken from
        frameSize  % Size of a single frame (2x1 vector: depth x AP length)
        scanArea   % Area that is scanned (depth x AP length)
        scanRes    % Scanning resolution (depth x AP length)
        movieLen   % Duration of movie in # of frames
        dt         % Movie Sampling rate
        stats      % Place to store basic statistics about the fu movie
%         meta   = struct('filename',[],'dataPath',[],'savePath',[],'subjectID',[]);
        errs   = struct('burst', [], 'motion', []);       % Error statistics
        post   = struct('movDenoised', []);                % post-processing results
        PCAout = struct();                                % PCA output
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Below are the set of methods associated with this class

    methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Main function
        function uo = fuobj(movieSpec, varargin)
            p = inputParser;                                               % Set up an object to parse all of the various inputs
            p.addParameter('nFrames'    , []);                             % Can opt not to load the entire movie
            p.addParameter('scanArea'   , []);                             % Need to input either the scan area...
            p.addParameter('scanRes'    , []);                             % ...or scan resolution
            p.addParameter('dt'         , []);                             % Option to give the temporal scan resolution...
            p.addParameter('frameTimes' , []);                             % ...and/or the times of frames.
            p.addParameter('events'     , []);                             % Option to input an event vector. If specified as a non-empty array will override the loaded event vector
            parse(p,varargin{:});
            p = p.Results;
            
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
                uo  = checkFrameTimes(uo, p);                              % Check the parameters for frame rate and time
                uo            = makeDataAndSavePath(uo, []);               % ...Create the different directories needed for saving
            elseif isstring(movieSpec)||ischar(movieSpec)                  % If a string (path of the movie) was passed in...
                assert(isfile(movieSpec))                                  % ...Make sure the file actually exists
                assert(strcmp(movieSpec(end-3:end),'.mat'))                % ...For now only mat files are permitted
                uo = loadUltrasoundDataFromMat(uo, movieSpec, p.nFrames);  % ...Load the movie and events vector
                movSz        = size(uo.movie);                             % ...Get the size of the movie
                uo.frameSize = movSz(1:2);                                 % ...Get the size of each frame
                uo.movieLen  = movSz(3);                                   % ...Get the movie duration
                uo.eventVec  = checkEventVec(uo.eventVec, p.events);
                uo           = checkFrameTimes(uo, p);
%                 if ~isempty(p.events); uo.eventVec = p.events; end        % Pass through the event vector
                uo           = makeDataAndSavePath(uo, movieSpec);         % ...Extract directory and file information, and create
            % elseif isequal(class(movieSpec), 'function_handle')
            else
                error('Cannot recognize movieSpec: invalid entry!')
            end
            
            uo = setResolution(uo, p, uo.frameSize);
            
        end
        function uo = makeDataAndSavePath(uo, pathFull)
            if isempty(pathFull)
                uo.dataPath  = [];                                         % The current data path is unknown! 
                uo.fileName  = [];                                         % The current data file name is unknown! 
                saveTmp      = pwd;                                        % Create a potential save name that follows the basic convention
                saveTmp      = [saveTmp,'/',datestr(now,'yyyymmdd'),...
                                                             'fuResults']; % Make a subfolder based on the date
                uo.savePath  = saveTmp;                                    % save the subfolder to the object
                if ~isfolder(saveTmp)
                    mkdir(saveTmp)
                end
            else
                [uo.dataPath, uo.fileName] = fileparts(pathFull);          % Extract the data path and file name
                saveTmp = strrep(uo.dataPath,'Data','Analysis');           % Create a potential save name that follows the basic convention
                if isfolder(saveTmp)
                    uo.savePath = saveTmp;                                 % Set the directory which all the results are saved to
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
       varargout = getBasicStats(fuobj, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      This function finds burst errors (frames with excessively high
%      error rates):
       varargout = findBurstFrames(fuobj, varargin);                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      This function finds motion correction errors (frames where the
%      motion correction failed):
       varargout = findMotionCorrectionError(fuobj, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      This function plots examples of found burst errors. Must run
%      'findBurstFrames' above first (or provide a vector of logical
%      values for fuobj.
       plotBurstErrors(fuobj, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        displayMovie(fuObj, varargin)                                      % Display the movie using MovieSlider
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        displayExampleFrame(fuObj, varargin)                               % Plot an example frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        displayExampleTraces(fuObj, varargin)                              % Plot example time traces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        displayMotionError(fuObj, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        denoiseTraces(fuObj, varargin)                                     % Wavelet denoising of traces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fuPCA(fuObj, varargin)                                             % Run PCA on the data
        % EventShape = event2timeseries(fuObj, varargin)
        % [eventframe, idx] = event2frame(fuObj, events, varargin)
        eSTH = eventSTH(fuObj, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       varargout = computeCorrsWithBaseline(fuobj, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    end
end



function uo = loadUltrasoundDataFromMat(uo, fileName, varargin)

% mov = loadUltrasoundDataFromMat(fileName, varargin)
%
% Function to load up an ultrasound dataset
%
% 2019 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin > 2
    nFrames = varargin{1};
else
    nFrames = [];
end

if nargin > 3
    nameToLoad = varargin{2};
else
    nameToLoad = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data

m = matfile(fileName);                                                     % Create mat file object to load data from

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
    nameToLoad = fNames{max3D};                                            % Pull out the name of the biggest 3D object
end

if isempty(nFrames); nFrames = fieldSize{max3D,1}(3);  end                 % If nFrames is empty then load all the frames

if numel(nFrames) == 1
    uo.movie = m.(nameToLoad)(:,:,1:nFrames);
elseif numel(nFrames) == 2
    uo.movie = m.(nameToLoad)(:,:,min(nFrames):max(nFrames));
elseif numel(nFrames) > 2
    uo.movie = m.(nameToLoad)(:,:,nFrames);
else
    error('Invalid frame selection!')
end

if any(strcmp(fNames,'events'));    uo.eventVec = m.events;
else;                               uo.eventVec = [];
end

if any(strcmp(fNames,'frames'));    uo.frameTimes = m.frames;
else;                               uo.frameTimes = [];
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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