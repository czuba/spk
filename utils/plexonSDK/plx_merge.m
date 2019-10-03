%   Merges 2 PLX files.
%   Usage:
%   mi = plx_reader(filename1,filename2, [filenameOut])
%     Optional input for output file name;
%     Returns mergeindex for use with plx_split.
% 
%   e.g.
%   mi = plx_merge('file1.plx',file2.plx','outfile.plx');
%   creates a composite plx file named outfile.plx and
%   returns the mergeindex (mi) which is the split point
%   in units of 10E5 ticks
%   --if no output filename is provided, will create
%     default name as filename1+filename2.plx