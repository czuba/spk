function outPath = homeTilda(inPath)

if isunix % mac or linux
    [~,homeDir] = unix('echo $HOME');
    homeDir = homeDir(1:end-1);

    if contains(inPath, homeDir)
        i = strfind(inPath, homeDir);
        inPath(i:length(homeDir)) = [];
        outPath = fullfile('~',inPath);
    else
        % do nothing
        outPath = inPath;
    end
    
else
    % Windoze...do nothing
    outPath = inPath;
end

