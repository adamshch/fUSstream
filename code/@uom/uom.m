classdef uom < handle
    
    
%    This file defines the functional ultrasound object fumulti
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        dataPath        % Path where all results are saved
        fileNames       % Name of file where data exists
        batchNames      % Name of file where data exists
        nFiles          % The number of files accessed
        subjectID       % Name of subject the movie was taken from
        movieLens       % Duration of movie in # of frames
        fileMetaData    % meta data for all the files accessed
        uo = cell(1);           % the actual uoects
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Below are the set of methods associated with this class

    methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Main function
        function uom = uom(dataPointer, varargin)
            p = inputParser;                                               % Set up an object to parse all of the various inputs
            p.addParameter('batchNames', []);                              % Optional names for the different datasets
            parse(p,varargin{:});
            p = p.Results;
            
            if isempty(dataPointer);   dataPointer = pwd; end              % If nothing is given, treat the current folder as the main path.
            
            if iscell(dataPointer)                                         % If a cell array (i.e., the actual movies) are passed in...
                uom.dataPath  = pwd;
                uom.nFiles    = numel(dataPointer);
                uom.fileNames = cell(uom.nFiles,1);                        % No file names if the data is put in directly!
                uom           = loadMultiData(uom, dataPointer);           % Load the data objects directly
            elseif isstring(dataPointer)||ischar(dataPointer)              % If a string (path to the folder with the movies) was passed in...
                uom.dataPath     = dataPointer;
%                 uom.fileMetaData = dir([dataPointer, '**/*.mat']);         % Get all the MAT files in the directory and sub-directories
                uom.fileMetaData = pickAllFusFiles(dataPointer);

                uom.nFiles       = numel(uom.fileMetaData);
                uom              = loadMultiMat(uom, uom.fileMetaData);
            else
                error('Cannot recognize movieSpec: invalid entry!')
            end
            uom = nameBatches(uom,p.batchNames);                           % Set up names of datasets
            uom.calcAllBasicStats();
            uom.removeTrivialDatasets(); 
        end                                                                % END OF MAIN FUNCTION
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Additional methods in other functions go here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Motion correction functions
        uom = multiMotionTest(uom, varargin); 
        uom = multiMotionCorrect(uom, varargin)
        
        uom = multiEventCorr(uom, varargin)
        uom = ensureAllMasks(uom, varargin)
        uom = calcAllBasicStats(uom)
        uom = removeTrivialDatasets(uom)
    end                                                                    % END OF METHODS
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function fileData = pickAllFusFiles(dataPointer)

fileData     = dir([dataPointer, '**/*.mat']);                             % Get all the MAT files in the directory and sub-directories
% motionFiles  = dir([dataPointer, '**/*motionCorrected.mat']);              % Get all the MAT files in the directory and sub-directories thhat are motion corrected files
% denoiseFiles = dir([dataPointer, '**/*denoised.mat']);                     % Get all the MAT files in the directory and sub-directories that are denoised data files

% fileData = setdiff(fileData, motionFiles);                                 % Remove all the motion files
% fileData = setdiff(fileData, denoiseFiles);                                % Remove all the denoised files

pickFiles = true(size(fileData));

for ll = 1:numel(fileData)
    if numel(fileData(ll).name) > 19
        if strcmp('motionCorrected.mat', fileData(ll).name(end-18:end))    % Remove all the motion files
            pickFiles(ll) = false;
        elseif strcmp('denoised.mat', fileData(ll).name(end-11:end))       % Remove all the denoised files
            pickFiles(ll) = false;
        elseif strcmp('mask.mat', fileData(ll).name(end-7:end))            % Remove all the mask files
            pickFiles(ll) = false;
        end
    end
    if pickFiles(ll)
        m         = matfile([fileData(ll).folder,'/',fileData(ll).name]);
        fNames    = fieldnames(m);                                         % Get all the field names
        fieldSize = cell(numel(fNames),2);                                 % Initialize a cell array to store sizes of all mat variables
        size3d    = zeros(numel(fNames),1);
        for kk = 1:numel(fNames)
            fieldSize{kk,1} = size(m.(fNames{kk}));                        % Get the size of all objects in the mat file
            fieldSize{kk,2} = numel(fieldSize{kk,1});                      % Get the dimensions of all objects in the mat file
            if fieldSize{kk,2}==3
                size3d(kk) = fieldSize{kk,1}(3);                           % Get the dimensions of all objects in the mat file
            else 
                size3d(kk) = 0;                                            % Get the dimensions of all objects in the mat file
            end
        end
        idx3D  = cell2mat(fieldSize(:,2)) == 3;                            % Find which arrays have 3 dimensions
        idx3D  = idx3D&(size3d>100);                                       % Find which arrays also have more than 100 frames
        if ~any(idx3D)
            pickFiles(ll) = false;
        end
    end
end

fileData = fileData(pickFiles);

end

function uom = loadMultiData(uom, fileData)

% mov = loadUltrasoundDataFromMat(fileName, varargin)
%
% Function to load up an ultrasound dataset
%
% 2019 - Adam Charles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data

for ll = 1:numel(fileData)
    uom.uo{ll} =  uo(fileData{ll}, 'nFrames', []);                   % Load the dataset
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function uom = loadMultiMat(uom, fileData)

% mov = loadUltrasoundDataFromMat(fileName, varargin)
%
% Function to load up an ultrasound dataset
%
% 2019 - Adam Charles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
uoCell = cell(1,numel(fileData));
for ll = 1:numel(fileData)
    uoCell{ll} =  uo([fileData(ll).folder,'/',fileData(ll).name], ...
                                                          'nFrames', []);  % Load the dataset
end
uom.uo = uoCell;
clear uo

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to name batches

function uom = nameBatches(uom,batchNames)

if isempty(batchNames)
    uom.batchNames = cell(uom.nFiles,1);
    for ll = 1:uom.nFiles
        uom.batchNames = sprintf('batch%2d',ll);
    end
elseif iscell(batchNames)&&(numel(batchNames) == uom.nFiles)
else
    warning('Bad input for batchNames: using default naming')
    uom = nameBatches(uom,[]);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%