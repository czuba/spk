%% WF plotting example:

spk = dvK.uprb;

figure;
nUnits = size(spk.wf.mu,3);
spx = ceil(sqrt(nUnits));
showCi = 0;
yl = min(range(spk.wf.mu, 'all'), 800)/2*[-1 1];

for u = 1:nUnits
    % Plot it
    h = subplot(spx, spx, u);
    plot(h, spk.wf.mu(:,:,u)); 
    if showCi
        hold on
        jnk = spk.wf.ci(:,:,:,u);
        plot(h, jnk(:,:), '--', 'yliminclude','off');
    end
    box off
    ylim(yl)
        
    % Label it
    if isfield(spk, 'sortId')
        titl = sprintf('\n%s',spk.sortId{u,2});
    else
        titl = '';
    end
    title(sprintf('unit %03d (%d)%s',u, length(spk.n(u)), titl), 'fontsize',8)

end
