function unstickPopup
% function unstickPopup
%
% Scrappy code snippet to unstick a Matlab popup window (function help blurbs are most frequent culprets...2019a)
%
% Code from:  https://www.mathworks.com/matlabcentral/answers/514677-how-can-i-clear-stuck-help-popups
% 
% 2020-02-01  TBC  Pasted into function.
% 

ws = com.mathworks.mwswing.window.MJFullWindowRegistry.windows;
while ws.hasMoreElements
    w=ws.nextElement;
    if any(strcmp(w.getName,{'HelpPopup','FunctionHints'}))
        w.dispose;
    end
end

end % main function