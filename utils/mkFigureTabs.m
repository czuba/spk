function figH = mkFigureTabs(nTabs, grpName)

if nargin<2 || isempty(grpName)
    grpName = 'tabGroup';
end

if isscalar(nTabs)
    tabNums = 1:nTabs;
else
    tabNums = nTabs;
end
nTabs = length(tabNums);

desktop = com.mathworks.mde.desk.MLDesktop.getInstance; %#ok<JAPIMATHWORKS>
myGroup = desktop.addGroup(grpName);
desktop.setGroupDocked(grpName, 0);

myDim   = java.awt.Dimension(1, nTabs);
% 1: Maximized, 2: Tiled, 3: Floating
desktop.setDocumentArrangement(grpName, 1, myDim)

figH    = gobjects(1, nTabs);

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
warning('off', 'MATLAB:ui:javaframe:PropertyToBeRemoved');

for iFig = 1:nTabs
   figH(iFig) = figure(tabNums(iFig));  %'WindowStyle', 'docked');%, ...
      %'Name', sprintf('Figure %d', iFig));%, 'NumberTitle', 'off');
      clf;
      set(figH(iFig), 'WindowStyle', 'docked');
%    drawnow;
%    pause(0.02);  % Magic, reduces rendering errors
   set(get(handle(figH(iFig)), 'javaframe'), 'GroupName', grpName); %#ok<*JAVFM>
%    plot(1:10, rand(1, 10));
end
