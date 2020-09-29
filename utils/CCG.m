function [ccg,t,tau,c] = CCG(times,id,varargin)

%CCG - Compute multiple cross- and auto-correlograms, or cross-covariances
%
%  USAGE
%
%    [ccg,t,tau,c] = CCG(times,id,<options>)
%
%    times          times of all events (see NOTE below)
%    id             ID for each event (e.g. unit ID) from 1 to n
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'binSize'     bin size in s (default = 0.01)
%     'duration'    duration in s of each xcorrelogram (default = 2)
%     'nBins'       number of bins (default = duration/binSize)
%     'smooth'      smoothing size in bins (0 = no smoothing, default)
%     'groups'      group number (1 or 2) for each event, used to restrict
%                   cross-correlograms to pairs across two groups of events
%                   (see EXAMPLE #2 below)
%     'mode'        'ccg' or 'ccv' (default = 'ccg')
%     'alpha'       significance level to determine correlated pairs
%     'totalTime'   total recording duration in s (if different from the
%                   default = max(times) - min(times))
%    =========================================================================
%
%  NOTES
%
%    The size of the cross-correlograms can be supplied either as a duration
%    (property 'duration') or as an number of bins (property 'nBins').
%
%  OUTPUT
%      ccg          value of cross-correlograms or cross-covariances
%                   dimensions are (nbins,m,n) where m is the number of
%                   reference time series (e.g. reference units) and n the
%                   number of referenced time series (in general m = n,
%                   except when using option 'groups')
%      t            time bins
%      tau          lag times for a each pair (mode 'ccv' only)
%      c            maximum cross-covariance for a each pair (mode 'ccv' only)
%
%
%  NOTE
%
%    Parameters 'times', 'id' and 'group' can be obtained using <a href="matlab:help CCGParameters">CCGParameters</a>.
%    As a special case, when computing the correlograms of spike trains, one
%    can use the output of <a href="matlab:help GetSpikes">GetSpikes</a> either directly or in combination with
%    <a href="matlab:help CCGParameters">CCGParameters</a>. See EXAMPLES below.
%
%  EXAMPLES
%
%    % Auto- and cross-correlograms between all neurons
%    spikes = GetSpikes('output','numbered');
%    [ccg,t] = CCG(spikes(:,1),spikes(:,2));
%
%    % Only tetrode #1 vs tetrode #2 (e.g. mPFC vs HPC neurons)
%    pfc = GetSpikes([1 -1],'output','numbered');
%    hpc = GetSpikes([2 -1],'output','numbered');
%    [s,ids,groups] = CCGParameters(pfc,hpc,2);
%    [ccg,t] = CCG(s,ids,'groups',groups);
%
%    % Between stimulations and MUA spikes
%    spikes = GetSpikes;
%    stimulatios = GetEvents('Stimulation');
%    d = [spikes(:,1) ones(size(spikes,1)) ; stimulations 2*ones(size(stimulations,1))];
%    d = sortrows(d);
%    [ccg,t] = CCG(d(:,1),d(:,2));
%
%    % To compute cross-covariances
%    [ccv,t,tau,C] = CCG(times,ids,'mode','ccv');
%
%  SEE
%
%    See also CCGParameters, ShortTimeCCG.

% Copyright (C) 2012-2013 by Michaël Zugaro, Marie Goutierre
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Default values
d = 2;
duration = [];
binSize = 0.01;
nBins = [];
smooth = 0;
groups = [];
mode = 'ccg';
alpha = 0.05;

% Check parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdvector(times),
	error('Parameter ''times'' is not a real-valued vector (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdscalar(id) && ~isdvector(id),
	error('Parameter ''id'' is not a real-valued scalar or vector (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdscalar(id) && length(times) ~= length(id),
	error('Parameters ''times'' and ''id'' have different lengths (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
id = id(:);
times = times(:);
totalTime = max(times)-min(times);

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help CCG">CCG</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'binsize',
			binSize = varargin{i+1};
			if ~isdscalar(binSize,'>0'),
				error('Incorrect value for property ''binSize'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'nbins',
			nBins = varargin{i+1};
			if ~isiscalar(nBins,'>0'),
				error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'duration',
			duration = varargin{i+1};
			if ~isdscalar(duration,'>0'),
				error('Incorrect value for property ''duration'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'smooth',
			smooth = varargin{i+1};
			if ~isdscalar(smooth,'>0'),
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'groups',
			groups = varargin{i+1};
			if ~isempty(groups) && ~isdvector(groups) && length(times) ~= length(groups)
				error('Incorrect value for property ''groups'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'alpha',
			alpha = varargin{i+1};
			if ~isdscalar(alpha,'>0'),
				error('Incorrect value for property ''alpha'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'mode',
			mode = varargin{i+1};
			if ~isastring(mode,'ccg','ccv'),
				error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'totaltime',
			totalTime = varargin{i+1};
			if ~isdscalar(totalTime,'>0'),
				error('Incorrect value for property ''totaltime'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).']);
	end
end

% Determine binSize/duration
if isempty(nBins),
	if isempty(duration),
		duration = d;
	end
else
	if isempty(duration),
		duration = nBins*binSize;
	elseif duration ~= binSize*nBins,
		error('Incompatible ''duration'' and ''nBins'' parameters (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
	end
end

tau = [];
c = [];

% Number of IDs, number of bins, etc.
if length(id) == 1,
	id = ones(length(times),1);
	nIDs = 1;
else
	nIDs = length(unique(id));
end
if nIDs ~= max(id),
	error('Incorrect IDs (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
halfBins = round(duration/binSize/2);
nBins = 2*halfBins+1;
t = (-halfBins:halfBins)'*binSize;

if length(times) <= 1,
	return
end

% Sort events in time and compute CCGs
[times,i] = sort(times);
id = id(i);
if ~isempty(groups),
	groups = groups(i);
end
counts = CCGEngine(times,id,binSize,halfBins);

% Reshape the results
n = max(id);
counts = reshape(counts,[nBins n n]);
if n < nIDs,
	counts(nBins,nIDs,nIDs) = 0;
end

% Restrict the results to inter-group CCGs if requested
if ~isempty(groups),
	group1 = unique(id(groups == 1));
	group2 = unique(id(groups == 2));
	nGroup1 = length(group1);
	nGroup2 = length(group2);
	ccg = zeros(nBins,nGroup1,nGroup2);
	for i = 1:nGroup1,
		for j = 1:nGroup2,
			ccg(:,i,j) = Smooth(flipud(counts(:,group1(i),group2(j))),smooth);
		end
	end
else
	ccg = zeros(nBins,nIDs,nIDs);
	% Compute corr(A,B) for each unique unordered pair (A,B)
	for g1 = 1:nIDs,
		for g2 = g1:nIDs
            if all(~smooth)
                ccg(:,g1,g2) = flipud(counts(:,g1,g2));
            else
                ccg(:,g1,g2) = Smooth(flipud(counts(:,g1,g2)),smooth);
            end
		end
	end
	% corr(B,A) and corr(B,A) symmetric
	for g1 = 1:nIDs,
		for g2 = 1:g1-1,
			ccg(:,g1,g2) = flipud(squeeze(ccg(:,g2,g1)));
		end
	end
end


if strcmp(mode,'ccv'),

	% Determine mean event rate for each ID
	eventRate = zeros(nIDs,1);
	for i = 1:nIDs,
		eventRate(i) = sum(id==i)/totalTime;
	end

	% Determine standardized cross-covariances
	ccv = zeros(size(ccg));
	tau = zeros(size(ccg,2),size(ccg,3));
	c = zeros(size(ccg,2),size(ccg,3));

	nPairs = size(ccg,2)*size(ccg,3);
	disp(['# pairs: ' int2str(nPairs)]);

	threshold = sqrt(2)*erfinv(1-(alpha/length(t)));

	for i = 1:size(ccg,2),
		for j = 1:size(ccg,3),
		
			% Compute and normalize CCVs from CCGs
			if ~isempty(groups),
				rate = eventRate(group1(i))*eventRate(group2(j));
			else
				rate = eventRate(i)*eventRate(j);
			end
			ccv(:,i,j) = sqrt((binSize*totalTime)/rate) * (ccg(:,i,j)/(binSize*totalTime)-rate);

			% Smooth with a 3-bin boxcar
			data = ccv(:,i,j);
			top = flipud(data(1:size(ccg,1),:));
			bottom = flipud(data(end-size(ccg,1)+1:end,:));
			data = [top;data;bottom];
			data = filter([1 1 1],3,data);
			n = size(data,1);
			d = n - size(ccg,1);
			start = d/2+1;
			stop = start + size(ccg,1) - 1;
			ccv(:,i,j) = data(start:stop);

			% Find the peak lag time and value
			[~,maxIndex] = max(ccv(:,i,j));
			tau(i,j) = median(t(maxIndex));
			c(i,j) = median(ccv(maxIndex,i,j));
			% Previous version of the code (discard?)
			% [~,maxIndex] = max(ccv(:,i,j));
			% tau(i,j) = t(maxIndex);
			% c(i,j) = median(ccv(max(1,maxIndex-3):min(end,maxIndex+3),i,j));

			% Keep only significantly correlated pairs
			if ~any(abs(ccv(:,i,j))>threshold),
				tau(i,j) = NaN;
			end
			
		end
	end

	nCorrelatedPairs = sum(~isnan(tau(:)));
	disp(['# significantly correlated pairs: ' int2str(nCorrelatedPairs)]);

	ccg = ccv;
	
end


end %main function


%% Sub-Functions
%
% FMAToolbox dependencies
% 

%% isdvector.m
function test = isdvector(x,varargin)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help isdvector">isdvector</a>'' for details).');
end

% Test: double, vector
test = isa(x,'double') & isvector(x);

% Ignore NaNs
x = x(~isnan(x));

% Optional tests
for i = 1:length(varargin),
	try
		if varargin{i}(1) == '#',
			if length(x) ~= str2num(varargin{i}(2:end)), test = false; return; end
		elseif isstring(varargin{i},'>','>=','<','<='),
			dx = diff(x);
			if ~eval(['all(0' varargin{i} 'dx);']), test = false; return; end
		else
			if ~eval(['all(x' varargin{i} ');']), test = false; return; end
		end
	catch err
		error(['Incorrect test ''' varargin{i} ''' (type ''help <a href="matlab:help isdvector">isdvector</a>'' for details).']);
	end
end

end %isdvector


%% isdscalar.m
function test = isdscalar(x,varargin)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help isdscalar">isdscalar</a>'' for details).');
end

% Test: double, scalar
test = isa(x,'double') & isscalar(x);

% Ignore NaN
x = x(~isnan(x));

% Optional tests
for i = 1:length(varargin),
	try
		if ~eval(['x' varargin{i} ';']), test = false; end
	catch err
		error(['Incorrect test ''' varargin{i} ''' (type ''help <a href="matlab:help isdscalar">isdscalar</a>'' for details).']);
	end
end

end %isdscalar


%% Smooth.m
function smoothed = Smooth(data,smooth)

%  vector = min(size(data)) == 1;
vector = isvector(data);
matrix = (~vector & length(size(data)) == 2);
if ~vector & ~matrix,
	error('Smoothing applies only to vectors or matrices (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).');
end

% Vectors must be 'vertical'
if size(data,1) == 1,
	data = data';
end

if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).');
end

% if ~isdvector(smooth,'>0') | (matrix & length(smooth) > 2) | (vector & length(smooth) ~= 1),
% 	error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).');
% end

% If Sh = Sv = 0, no smoothing required
if all(smooth==0),
	smoothed = data;
	return
end

if length(smooth) == 1,
	% For 2D data, providing only one value S for the std is interpreted as Sh = Sv = S
	smooth = [smooth smooth];
end

% Build Gaussian kernels
[vSize,hSize] = size(data);
% 1) Vertical kernel
vKernelSize = min([vSize 1001]);
r = (-vKernelSize:vKernelSize)'/vKernelSize;
vKernelStdev = smooth(1)/vKernelSize;
vKernel = exp(-r.^2/(vKernelStdev+eps)^2/2);
vKernel = vKernel/sum(vKernel);
% 2) Horizontal kernel
hKernelSize = min([hSize 1001]);
r = (-hKernelSize:hKernelSize)/hKernelSize;
hKernelStdev = smooth(2)/hKernelSize;
hKernel = exp(-r.^2/(hKernelStdev+eps)^2/2);
hKernel = hKernel/sum(hKernel);
if vector,
	% Vector smoothing
	% Prepend/append data to limit edge effects
	top = flipud(data(1:vKernelSize));
	bottom = flipud(data(end-vKernelSize+1:end));
	data = [top;data;bottom];
	% Convolve (and return central part)
	tmp = conv(vKernel,data);
	n = size(tmp,1);
	d = n - vSize;
	start = d/2+1;
	stop = start + vSize - 1;
	smoothed = tmp(start:stop,:);
else
	% Matrix smoothing
	% Convolve
	if smooth(1) == 0,
		% Smooth only across columns (Sv = 0)
		% Prepend/append data to limit edge effects
		left = fliplr(data(:,1:hKernelSize));
		right = fliplr(data(:,end-hKernelSize+1:end));
		data = [left data right];
		for i = 1:size(data,1),
			tmp = conv(hKernel,data(i,:));
			n = size(tmp,2);
			d = n - hSize;
			start = d/2+1;
			stop = start + hSize - 1;
			smoothed(i,:) = tmp(:,start:stop);
		end
	elseif smooth(2) == 0,
		% Smooth only across lines (Sh = 0)
		% Prepend/append data to limit edge effects
		top = flipud(data(1:vKernelSize,:));
		bottom = flipud(data(end-vKernelSize+1:end,:));
		data = [top;data;bottom];
		for i = 1:size(data,2),
			tmp = conv(vKernel,data(:,i));
			n = size(tmp,1);
			d = n - vSize;
			start = d/2+1;
			stop = start + vSize - 1;
			smoothed(:,i) = tmp(start:stop);
		end
	else
		% Smooth in 2D
		% Prepend/append data to limit edge effects
		top = flipud(data(1:vKernelSize,:));
		bottom = flipud(data(end-vKernelSize+1:end,:));
		data = [top;data;bottom];
		left = fliplr(data(:,1:hKernelSize));
		right = fliplr(data(:,end-hKernelSize+1:end));
		data = [left data right];
		tmp = conv2(vKernel,hKernel,data,'same');
		n = size(tmp,1);
		d = n - vSize;
		vStart = d/2+1;
		vStop = vStart + vSize - 1;
		n = size(tmp,2);
		d = n - hSize;
		hStart = d/2+1;
		hStop = hStart + hSize - 1;
		smoothed = tmp(vStart:vStop,hStart:hStop);
	end
end

end


