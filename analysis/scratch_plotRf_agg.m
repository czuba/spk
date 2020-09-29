function H = scratch_plotRf_agg(dvRf, saveout, withDate)
% function H = scratch_plotRf_agg(dvRf, saveout, withDate)
% 
% 

if nargin<3 || isempty(withDate)
    withDate = 1;
end
if nargin<2 || isempty(saveout)
    saveout = 0;
end

nunits = size(dvRf.rf.fit,1);

% site name
pdsName = dvRf.info.pdsName;
if isfield(dvRf.paths, 'site')
    [~,siteName] = fileparts(dvRf.paths.site);
else
    siteName = pdsName;
end


% xy grid for plotting fits
spc = min([diff(dvRf.rf.xs); diff(dvRf.rf.ys)]);
[xxRF, yyRF] = meshgrid(dvRf.rf.xs(1)-10:.2:dvRf.rf.xs(end)+10, dvRf.rf.ys(1)-10:.2:dvRf.rf.ys(end)+10);

rfx = dvRf.rf.fit(:,3);
rfy = dvRf.rf.fit(:,4);
rfw = dvRf.rf.fit(:,5)*2.355; % fwhm from sd

% default stereo probe depths
trodeDepth = fliplr(kron(0:100:1500, [1,1]))';

% try to get depth from file name (may fail)
fnParts = strsplit(dvRf.info.pdsName, {'_','-'});
if ~contains(fnParts{end-1}, {'M','P','L','A'})
    % clean it up
    dstr = fnParts{end-1};
    ii = [strfind(dstr,'h'), strfind(dstr,'m'), strfind(dstr,'-')];
    dstr(ii) = [];
    
    siteDepth = str2num(dstr);
else
    siteDepth = trodeDepth(1);
end

% depth of each trode recording site (in mm)
rfDepth = (siteDepth - trodeDepth)/1000;

% colormaps
% cm = colorcet('D11', 'N',nunits/2); % isoluminant colormap %colormap(winter)
cm = colorcet('R2', 'N',nunits/2); % rainbow colormap %colormap(winter)
cm = kron(cm,[1 1]');



% Plot rfFit contour
H = figureFS([],'portrait',[8,8]);
if isfield(dvRf.paths,'figs')
    figDir = dvRf.paths.figs;
else
    figDir = fullfile(dvRf.paths.root,'figs');
end
set(H, 'name', sprintf('%s_rfAgg',siteName), 'tag',figDir); %dvRf.paths.figs);

% add stim location grid to last image
plot( dvRf.rf.xyz(:,1), dvRf.rf.xyz(:,2), 'o','color',.32*[1 1 1], 'markersize',3)
hold on, axis equal
colormap(cm)
% c50 = .75*[1 1 1];
lw = .5;
fsz = 12;
dotMin = 20;

% rf centers
mkSz = max([dvRf.rf.r2; 0.4])
mkSz = (dvRf.rf.r2)./mkSz .*250+dotMin;
hp = scatter(rfx, rfy, mkSz, cm);
hp.MarkerEdgeColor = 'k';

% color depth
cl = [min(rfDepth), max(rfDepth)];
set(gca, 'clim', cl);
cb = colorbar('south');

xyunity
box off

titext = siteName;
if ~isequal(titext, pdsName)
    titext = {titext, pdsName};
end
    
title(titext, 'interpreter','none', 'fontsize',fsz)

% Axis limits
%   Same as single unit rf plots
% mapping grid lims
gl = [min(dvRf.rf.xyz); max(dvRf.rf.xyz)];
% make image pixels fall in corrrect locations
gl = gl + kron(diff(gl)./dvRf.rf.rate.condDims(1:3) ./2, [-1;1]); 

% plotting axes limits
% aa = axis .*1.2;
rfStimLims = gl(1:4) + 3*[-1,1,-1,1]; % make stimulus edges visible
% rfStimLims = gl(1:4) + 1*[-1,1,-1,1]; % make stimulus edges visible
% rfAxLims = [-3,22,-20,5];  % [-8,32,-25,10]; % a standard large area of screen
rfAxLims = [-4,24,-23,5];  
%     % Set axes lims inclusive of stim & standard
%     ii = [1,3];     rfAxLims(ii) = min([rfStimLims(ii); rfAxLims(ii)]);
%     ii = [2,4];     rfAxLims(ii) = max([rfStimLims(ii); rfAxLims(ii)]);
% 
% rfAxLims = rfStimLims;

tk = -35:5:30;
tkl = [""    "-30"    ""    "-20"    ""    "-10"    ""    "0"    ""    "10"    ""    "20"    ""    "30"];




set(gca, 'fontsize',fsz*1.5)

% Unit numbers
try
    % only if unsorted (too much variance on sorted to be helpful)
    if nunits==dv.uprb.info.rawInfo.nChannels
        nL = [1,8:8:32];
        
        noff = 0.6;
        % fsz = 18;
        for u = nL%1:nunits
            plot( rfx(u)*[1 1], rfy(u)*[1 1]+[0,noff], '-k', 'linewidth',1.32);
            ht = text(rfx(u), rfy(u)+noff+.1, num2str(u));    %sprintf('%2.1f',rfDepth(u)));
            set(ht, 'fontweight','bold', 'horizontalalignment','center',  'fontsize',fsz*1.5, 'verticalalignment','baseline', 'color',.02*[1 1 1])
        end
    end
end

% axis(aa)
axis(rfAxLims)
set(gca, 'XTick',tk, 'YTick',tk, 'XTickLabel',tkl, 'YTickLabel',tkl)
drawnow

if saveout
    set(H, 'name', sprintf('%s_rfAgg-1',siteName))
    saveFigTriplet(withDate, get(gcf,'name'), [0 1 0])
end

%% RF circles
fitcut = .01;%prctile(dvRf.rf.r2, 10);
brtLvl = .7; % brighten ammount
hc = [];

for u = 1:nunits
    if dvRf.rf.r2(u)<=fitcut
        ls_fit = '--';
    else
        ls_fit = '-';
    end
    c50 = brighten(cm(u,:),brtLvl);
    
    if isfield(dvRf.rf, 'fxn')
        fout = dvRf.rf.fxn(dvRf.rf.fit(u,:), {xxRF(:),yyRF(:)});
        rf50 = dvRf.rf.fit(u,1) + dvRf.rf.fit(u,2)/2;
        [~,hc(u)] = contour(xxRF, yyRF, reshape(fout, size(xxRF)), [1,1]*rf50, 'color',c50,'linestyle',ls_fit, 'linewidth',lw, 'xliminclude','off', 'yliminclude','off');
    end
end
uistack(hc, 'bottom')

if saveout
    set(H, 'name', sprintf('%s_rfAgg-2',siteName))
    saveFigTriplet(withDate, get(gcf,'name'), [0 1 0])
end

%% direction
% No...garbage mess of a plot!
% % %         fit2 = squeeze(dvXz.tune.fit2d.fit);
% % %         fit2r = dvXz.tune.fit2d.r2;
% % %         dirV = dvXz.tune.dir2d;
% % %         qlen = 5;
% % % 
% % %         hq = quiver(rfx', rfy', cos(fit2(4,:)).*fit2r.*qlen, sin(fit2(4,:)).*fit2r.*qlen, 0, 'filled')
% % %         % hq = quiver(rfx', rfy', cosd(dirV(1,:)).*dirV(2,:).*qlen, sind(dirV(1,:)).*dirV(2,:).*qlen, 0, 'filled')


