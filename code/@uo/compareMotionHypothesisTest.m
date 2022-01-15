function motTest = compareMotionHypothesisTest(uo, varargin)

% motTest = compareMotionHypothesisTest(uo, varargin)
%
% This function computes the motion comaprison metric (see the function 
% motionHypothesisTest()) both before and after motion correction and
% compares the results. All input arguments are passed directly into
% motionHypothesisTest().
%
% Optional parameters are:
%   'motionCorrect' - Select if the motion corrected data should be used: true/false (default false)                               
%   'denoised'      - Choose to display denoised movie true/false (default false) 
%   'plotOpt'       - Optional plotting flag. true/false (default false) 
%   'figNumber'     - Optional figure number to plot to (default 203) 
%   'motionCompare' - Choose which motion metric to correlate to (default 'total')
%   'useMask'       - Choose to use the mask and only consider "inside" pixels. true/false (default true)
%   'minHistNum'    - Choose to use the mask and only consider "inside" pixels (default 1000)
%   'distSelect'    - Choose the distance metric (default 'emd')
%   'shuffle'       - Choose whether to SHUFFLE the data to test the metric. true/false (default false)   
%
% Output is a struct with
%    motTest.preMotion -  motion metric per pixel before motion correction
%    motTest.postMotion - motion metric per pixel after motion correction
%    
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

motTest.preMotion  = uo.motionHypothesisTest('motionCorrect' , false, varargin{:});
motTest.postMotion = uo.motionHypothesisTest('motionCorrect' , true , varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First a histogram

edgeCut = 12;
distLims = 1.1*[0,max(vec(motTest.preMotion.distRatio))];

figure(1232)
if isinf(distLims(2))
    binLims = [linspace(distLims(1), 1e4, 99), Inf];
else
    binLims = linspace(distLims(1), distLims(2), 100);
end
if any(isnan(distLims))
    fprintf('Error with data: Nan values obtained for score values. Skipping plotting....\n')
else
    hPre = histcounts(vec(motTest.preMotion.distRatio(edgeCut:end-edgeCut+1, edgeCut:end-edgeCut+1)),  binLims);
    hPst = histcounts(vec(motTest.postMotion.distRatio(edgeCut:end-edgeCut+1, edgeCut:end-edgeCut+1)), binLims);

    subplot(2,2,[1,2]), bar(0.5*(binLims(1:end-1)+binLims(2:end)), [hPre(:),hPst(:)])

    box off; set(gca,'XLim',distLims, 'YScale','log', 'YLim', [0.5, 1e4])      % Remove box and set limits
    legend('Motion scores: no correction', 'Motion scores: with correction')   % Set the legend

    subplot(2,2,3), imagesc(motTest.preMotion.distRatio(edgeCut:end-edgeCut+1, edgeCut:end-edgeCut+1), distLims)
    pbaspect([uo.scanArea(2),uo.scanArea(1),1]);                               % Set the aspect ratio
    axis off;                                                              % Remove axes 
    title(sprintf('Score map: no correction'))                                 % Give a title to the image
    colorbar

    subplot(2,2,4), imagesc(motTest.postMotion.distRatio(edgeCut:end-edgeCut+1, edgeCut:end-edgeCut+1), distLims)
    pbaspect([uo.scanArea(2),uo.scanArea(1),1]);                               % Set the aspect ratio
    axis off;                                                                  % Remove axes 
    title(sprintf('Score map: with correction'))                               % Give a title to the image
    colorbar
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
