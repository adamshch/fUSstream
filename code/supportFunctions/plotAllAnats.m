function [uoMask, uoStatIms] = plotAllAnats(uom, p)

% [uoMask, uoStatIms] = plotAllAnats(uom, p)
% 
% Function to plot the anatomical images for sessions across multiple days
% 
% 2021 - Adam Charles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uoStatIms = cell(numel(uom.uo),3);

p.motionCorrect = 'true';
for ll = 1:numel(uom.uo)
    movieOut = uom.uo{ll}.makeGetMovieFunction(p);
    uoStatIms{ll,1}   = mean(movieOut,3);
    uoStatIms{ll,2}   = median(movieOut,3);
    uoStatIms{ll,3}   = prctile(movieOut,95,3);
end

[uoMask, uoStatIms] = getMotSizedMaskAndStats(uom, uoStatIms);

for ll = 1:numel(uom.uo)
    for kk = 1:size(uoStatIms,2)
        uoStatIms{ll,kk}(isnan(uoStatIms{ll,kk})) = 0;
    end
end



titleCell = {'06/28/2018', '06/29/2018', '08/22/2018', '08/24/2018'};
imgTypes  = {'Mean', 'Median', '95th pctle'};
[D1,D2] = size(uoStatIms);

figure(1)
for ll = 1:D1
    for kk = 1:D2
        subplot(D1,D2,kk + D2*(ll-1)), imagesc(uoStatIms{ll,kk}./max(abs(uoStatIms{ll,kk}(:))), [-0, 0.5])
        pbaspect([uom.uo{ll}.scanArea(2),uom.uo{ll}.scanArea(1),1]);   
        axis off; colormap gray 
        title(['Map:', imgTypes{kk},' ',titleCell{ll}])
    end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra functions

function [uoMask, uoStatIms] = getMotSizedMaskAndStats(uom, uoStatIms)


padSize = zeros(numel(uom.uo), 4);
for ll = 1:numel(uom.uo)
    padSize(ll,:) = uom.uo{ll}.post.motionPad;
end

padGlobal = max(padSize,[],1);
uoMask    = cell(1, numel(uom.uo));
for ll = 1:numel(uom.uo)
    padDiff = padGlobal - padSize(ll,:);
    uoMask{ll}  = padarray(padarray(uom.uo{ll}.mask, padGlobal([1,3]),0,'pre'), ...
                                              padGlobal([2,4]), 0, 'post');
    for kk = 1:size(uoStatIms,2)
        uoStatIms{ll,kk} = padarray(padarray(uoStatIms{ll,kk}, ...
                     padDiff([1,3]),0,'pre'), padDiff([2,4]), 0, 'post');
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
