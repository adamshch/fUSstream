%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Demo for using the Ultrasound Object class  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This demo contains a series of steps that show how to use the uo and uom  
% object classes. These classes are designed to make handling one or more
% ultrasound videos more managable by storing either full videos or links
% to saved files on disk. 
% 
% 2021 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial setup
% First the paths need to be added. Specifically all paths in 
% PATH-TO-US/code should be set. Either run from that directory or add the
% paths manually by changing the folder:

addpath(genpath('.'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Next get pointers to the actual data
% It's sometimes helpful when dealing with multiple datasets to have a path
% variable and a cell array referencing all the different datasets. This 
% block provides such an example.

globalObj.path     = '/home/adam/BigStorage/Data/fUs/';                    % Set path to the data
globalObj.dataSets = {'H134_LongPoisson_fUSi_062818_153132_fus.mat',...
                      'H134_LongPoisson_fUSi_082218_162003_fus.mat',...
                      'H134_LongPoisson_fUSi_092618_155625_pos1.mat', ...
                                       'DataExample_Codes/H134_data.mat'}; % The list of currently accessible files (will be updated)
globalObj.dataSelect = 4;                                                  % Choose which data to load (#2 is good for motion correction testing)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load a single dataset into a uo object
% uo requires the path and filename of the data bring loaded. 'nFrames' can
% subselect which and how many frames to load (empty means all frames). 
% 'loadToRAM' tells uo to load the enire movie rather than setting a 
% pointer to the data as a function that loades frames as needed. This is 
% better for single and small datasets but depends on the available
% computational resources.

fuo = uo([globalObj.path, globalObj.dataSets{globalObj.dataSelect}], ...
                     'nFrames', [],'preCalcStats',true,'loadToRAM',true);  % Load the dataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find burst errors
% One major effect of motion is the sudden burst of intensity across the
% entire image. This is caused by a small jerk on the probe. Identifying 
% these errors can be accomplished by removing frames in the tail of the
% frame intensity histogram and filling in the missing frames via linear
% interpolation. 

fuo.findBurstFrames()                                                      % Finds the burst error frames
fuo.findMotionCorrectionError();
[frameInPt,frameIDX] = fuo.computeBurstErrorInpainting();                  % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify motion correction errors
% Residual motion can be in the axial or lateral direction. Rigid motion
% can account for a significant portion of this and motion and can 
%

fuo.findMotionCorrectionError('batchSz', 5, 'useDenoised', false, ...
                           'searchBlock', 10, 'medBaseLine', 'firstmedian')% Example of motion correction on one dataset
fuo.correctResidualMotion('reCalc', true);                                 % One can force the moriotn correction to be recalculated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at the movie using MovieSlider
% This function requires the MovieSlider package available at 
%
%      https://github.com/sakoay/MovieSlider
%

fuo.displayMovie()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display some example traces from the movie. 

fuo.displayExampleTraces('numTraces', 10, 'pixSelect', 'rand', ...
         'normTraces', true, 'dispBurstErrs', true, 'dispMotErrs',true,...
                                                             'smoothLvl',6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display an example movie frame

fuo.displayExampleFrame('frameNo','rand','scaleFactor',[2,5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Denoise movie (in time) and view
% Wavelet denoising is applied trace-by-trace

fuo.denoiseTracesWavelet();                                                 % Denoise the data with default settings
fuo.displayMovie('denoised', 'true')                                        % Display the denoised movie using MovieSlider

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute PCA on the data

fuo.fuPCA('motionCorrect' , true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GraFT decomposition of a fUSi dataset

fuo.fuGraFT('n_dict', 45, 'lamForb', 0.3, 'lamCont', 0.4, 'lamCorr', 0.3, 'lamContStp', 0.9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Multiple ultrasound datasets
% To load multiple files, simply point the uom function to the folder that
% contains all the files to work with.

fuom = uom('/home/adam/BigStorage/Data/fUs/multiMotionTestFiles/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Motion correction across multiple files
% The uom object enables a single call to detect translational motion 
% across multiple sessions of the same plane. 

fuom.multiMotionTest('globalOrient', true, 'searchBlock' , 10, 'batchSz',...
                                              20, 'medBaseLine' , 'median')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
