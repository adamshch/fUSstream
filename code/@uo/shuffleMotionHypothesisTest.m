function motTest = shuffleMotionHypothesisTest(uo, varargin)

% motTest = shuffleMotionHypothesisTest(uo, varargin)
%
% This function 
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if isstr(varargin{1})
    motTest.noShuffle   = uo.motionHypothesisTest('motionCorrect', false, varargin{:}, 'shuffle', false);
    motTest.noShuffleMC = uo.motionHypothesisTest('motionCorrect', true,  varargin{:}, 'shuffle', false);
    motTest.withShuffle = uo.motionHypothesisTest('motionCorrect', false, varargin{:}, 'shuffle', true);
else
    motTest = varargin{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First a histogram

nBins   = 50;
edgeCut = 12;
mask2   = zeros(size(uo.mask));
mask2(edgeCut:end-edgeCut+1, edgeCut:end-edgeCut+1) = 1;
mask2 = mask2.*uo.mask;
distLims = 1.1*[0,max([max(vec(motTest.noShuffle.distRatio(mask2==1))), ...
    max(vec(motTest.noShuffleMC.distRatio(mask2==1))), ...
    max(vec(motTest.withShuffle.distRatio(mask2==1)))])];
mask2(mask2==0) = NaN;

figure(1233)
if isinf(distLims(2))
    binLims = [linspace(distLims(1), 1e4, nBins-1), Inf];
else
    binLims = linspace(distLims(1), distLims(2), nBins);
end
if any(isnan(distLims))
    fprintf('Error with data: Nan values obtained for score values. Skipping plotting....\n')
else
    hPre = histcounts(vec(motTest.noShuffle.distRatio(mask2==1)),   binLims);
    hMC  = histcounts(vec(motTest.noShuffleMC.distRatio(mask2==1)), binLims);
    hPst = histcounts(vec(motTest.withShuffle.distRatio(mask2==1)), binLims);
    
    subplot(2,3,[1,3]), bar(0.5*(binLims(1:end-1)+binLims(2:end)), [hPre(:),hMC(:),hPst(:)])

    box off; set(gca,'XLim',distLims, 'YScale','log', 'YLim', [0.5, 1e4])  % Remove box and set limits
    legend('Motion scores: Pre MC', 'Motion scores: Post MC', 'Motion scores: shuffle')          % Set the legend

    subplot(2,3,4), imagesc(mask2.*motTest.noShuffle.distRatio, 'AlphaData',~isnan(mask2), distLims)
    pbaspect([uo.scanArea(2),uo.scanArea(1),1]);                           % Set the aspect ratio
    axis off;                                                              % Remove axes 
    title(sprintf('Score map: no motion correction'))                                % Give a title to the image
    colorbar

    subplot(2,3,5), imagesc(mask2.*motTest.noShuffleMC.distRatio, 'AlphaData',~isnan(mask2), distLims)
    pbaspect([uo.scanArea(2),uo.scanArea(1),1]);                           % Set the aspect ratio
    axis off;                                                              % Remove axes 
    title(sprintf('Score map: motion corrected'))                                % Give a title to the image
    colorbar

    subplot(2,3,6), imagesc(mask2.*motTest.withShuffle.distRatio, 'AlphaData',~isnan(mask2), distLims)
    pbaspect([uo.scanArea(2),uo.scanArea(1),1]);                           % Set the aspect ratio
    axis off;                                                              % Remove axes 
    title(sprintf('Score map: shuffle'))                                   % Give a title to the image
    colorbar
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%