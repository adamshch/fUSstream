function plotShuffleTestHistograms(motAnalysis)


for ll = 1:numel(motAnalysis)
    
    motTest  = motAnalysis{ll};
    edgeCut  = 3;
    numBins  = 100;
    distLims = 1.1*[0,max([vec(motTest.noShuffle.distRatio); vec(motTest.noShuffleMC.distRatio)])];
    
    if isinf(distLims(2))
        binLims = [linspace(distLims(1), 1e4, 99), Inf];
    else
        binLims = linspace(distLims(1), distLims(2), numBins);
    end
    if any(isnan(distLims))
        fprintf('Error with data: Nan values obtained for score values. Skipping plotting....\n')
    else
        figure(1232)
        
        hPre = histcounts(vec(motTest.noShuffle.distRatio(edgeCut:end-edgeCut+1, edgeCut:end-edgeCut+1)),  binLims);
        hPst = histcounts(vec(motTest.noShuffleMC.distRatio(edgeCut:end-edgeCut+1, edgeCut:end-edgeCut+1)), binLims);
        hSfl = histcounts(vec(motTest.withShuffle.distRatio(edgeCut:end-edgeCut+1, edgeCut:end-edgeCut+1)), binLims);
        
        subplot(2,2,ll), bar(0.5*(binLims(1:end-1)+binLims(2:end)), [hPre(:),hPst(:),hSfl(:)])

        box off; set(gca,'XLim',distLims,'YScale','log','YLim',[0.5,1e5])  % Remove box and set limits
        legend('Motion scores: no correction', 'Motion scores: with correction')   % Set the legend

%         subplot(2,2,3), imagesc(motTest.noShuffle.distRatio(edgeCut:end-edgeCut+1, edgeCut:end-edgeCut+1), distLims)
%         pbaspect([fuo.scanArea(2),fuo.scanArea(1),1]);                     % Set the aspect ratio
%         axis off;                                                          % Remove axes 
%         title(sprintf('Score map: no correction'))                         % Give a title to the image
%         colorbar
% 
%         subplot(2,2,4), imagesc(motTest.noShuffleMC.distRatio(...
%                   edgeCut:end-edgeCut+1, edgeCut:end-edgeCut+1), distLims)
%         pbaspect([fuo.scanArea(2),fuo.scanArea(1),1]);                     % Set the aspect ratio
%         axis off;                                                          % Remove axes 
%         title(sprintf('Score map: with correction'))                       % Give a title to the image
%         colorbar

        %h = figure(1232); save2pdf(sprintf('/home/adam/Dropbox/MotionMetric_batchTMP.pdf', ll), h)
    end
end


end