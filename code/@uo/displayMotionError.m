function displayMotionError(uo, varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('baseSelect'    , 'median');                                % Select what criteria to calculate the baseline image as
p.addParameter('basePrctile'   , 25      );                                % For the situation that 'prctile' is chosen, which percentile should the baseline image be calculated using
p.addParameter('pixSelect'     , 'vary'  );                                % Select which pixels to select: either an array or 'vary' or 'rand'
p.addParameter('numTraces'     , []      );                                % Select hoe many traces to select
p.addParameter('normTraces'    , false   );                                % Choose if each trace is normalized to its own max or not
p.addParameter('dispBurstErrs' , false   );                                % Select underlay of burst error locations
p.addParameter('dispMotErrs'   , true   );                                 % Select underlay of motion correction errors
p.addParameter('smoothLvl'     , 0       );                                % Select denoising level
p.addParameter('figNo'         , 105     );                                % Select figure number

parse(p,varargin{:});
p = p.Results;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up the variables

[motDep, motLat, maxMot] = getMotionVectors(uo, p);                        % Check for motion errors if needed

tt = uo.dt*((1:numel(motLat))-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Actually perform plotting

figure(p.figNo); clf;
figure(p.figNo), hold on;

h = fill([tt(:); flipud(tt(:))]',[0*tt(:);...
                                flipud(motDep(:))]','r');
set(h,'facealpha',0.3)
legendNames{1} = 'Depth motion offset';
h = fill([tt(:); flipud(tt(:))]',[0*tt(:);...
                                flipud(motLat(:))]','b');
set(h,'facealpha',0.3)
legendNames{2} = 'AP motion offset';
totalMotErr = sqrt(motDep(:).^2+motLat(:).^2);
plot(tt,totalMotErr,'k')
legendNames{3} = 'Total motion offset';
ylabel('Motion offset (px)')
xlabel('Time (s)')
set(gca, 'XLim', [tt(1),tt(end)])
box off
legend(legendNames)
figure(p.figNo), hold off;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to check and extract motion errors

function [motDep, motLat, maxMot] = getMotionVectors(uo, p)

if p.dispMotErrs
    if ~isfield(uo.errs,'motion')                                          % Check to make sure that the motion correction error check was run
        warning('No motion correction error detection output found. Running findMotionCorrectionError with default settings...\n')
        uo.findMotionCorrectionError();                                    % Run the motion correction error detection
    end 
    motDep = uo.errs.motion.depth;                                         % Extract the depth motion offsets
    motLat = uo.errs.motion.ap;                                            % Extract the lateral motion offsets

    if numel(motDep)~=uo.movieLen                                          % If motDep is not the right size, then the batch version of motion correction must have been used...
        motDep = repelem(motDep(:),uo.errs.motion.params.batchSz);         % In that case replecate the elements to reflect the batch processing
        motDep = motDep(1:uo.movieLen);                                    % Truncate to accound fot the last block potentially being incomplete
    end

    if numel(motLat)~=uo.movieLen                                          % If motLat is not the right size, then the batch version of motion correction must have been used...
        motLat = repelem(motLat(:),uo.errs.motion.params.batchSz);         % In that case replecate the elements to reflect the batch processing
        motLat = motLat(1:uo.movieLen);                                    % Truncate to accound fot the last block potentially being incomplete
    end

    motDep(isnan(motDep)) = 0;
    motLat(isnan(motLat)) = 0;

    maxMot = max(max(abs(motDep)), max(abs(motLat)));
    motDep = 0.8*motDep/maxMot;
    motLat = 0.8*motLat/maxMot;
else
    motLat = [];
    motDep = [];
    maxMot = NaN;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%