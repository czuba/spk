function [allts allwave] = plx_reader(OpenedFileName,arg)
%   Reads in the data from PLX files.
%based on the SDK provided by Plexon
%   Usage:
%   [TimeStamps,WaveForms] = plx_reader(filename,arg)
%   e.g.
%   1) To read all Time Stamps, ignores the waveforms:
%      TimeStamps = plx_reader('file.plx',0);
%   2) To read all Time Stamps and waveforms:
%      [TimeStamps,WaveForms] = plx_reader('file.plx',0)
%   3) To reads all Time Stamps and waveforms, except those with sort-code = 0
%      TimeStamps = plx_reader('file.plx',1)
%   4) To reads all Time Stamps, except those with sort-code = 0
%      [TimeStamps,WaveForms] = plx_reader('file.plx',1)

if nargin <2
    arg = 0;
end
[OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreThresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(OpenedFileName);
[tscounts, wfcounts, evcounts, slowcounts] = plx_info(OpenedFileName,1);

if arg == 1
    tscounts(1,:) = 0;
end

ss = find(sum(tscounts));
if Trodalness == 4
    sel = ss(logical([1 diff(fix((ss-2)/4))]));
else
    sel = ss;
end
tscountsn = zeros(size(tscounts));
tscountsn(:,sel) = tscounts(:,sel);
%tscounts = tscountsn;
clear sel ss;

% need the event chanmap to make any sense of these
[~,evchans] = plx_event_chanmap(OpenedFileName);
evcounts = evcounts(evchans<=32);
evchans  = evchans( evchans<=32);

allts = zeros(sum([tscountsn(:); evcounts(:)]), 3);
if nargin <2
    allwave = zeros(sum([tscountsn(:); evcounts(:)]), NPW*Trodalness);
end

[xx,yy] = find(tscountsn);

%spikes
c = 0;
for i = 1:length(xx)
    allts(     c+(1:tscountsn(xx(i),yy(i))),1)  = yy(i)-1;
    allts(     c+(1:tscountsn(xx(i),yy(i))),2)  = xx(i)-1;
    if nargout == 1
        [nts allts(c+(1:tscountsn(xx(i),yy(i))),3)] = plx_ts(OpenedFileName, yy(i)-1 , xx(i)-1 );
    else
        for j = 1:Trodalness
            if (tscounts(xx(i),(fix((yy(i)-2)/4)*4+j+1))~=0)
                [nts npw allts(c+(1:tscountsn(xx(i),yy(i))),3), allwave(c+(1:tscountsn(xx(i),yy(i))),(1:NPW)+(j-1)*NPW)] = plx_waves(OpenedFileName, fix((yy(i)-2)/4)*4+j , xx(i)-1 );
            end
        end
    end
    c = c+nts;
end


%events

for iev = 1:length(evchans)
    if ( evcounts(iev) > 0 )
        evch = evchans(iev);
        [nevs, allts(c+(1:evcounts(iev)),3)] = plx_event_ts(OpenedFileName, evch);
        allts(c+(1:evcounts(iev)),1) = 0;
        allts(c+(1:evcounts(iev)),2) = evchans(iev);
        c = c + nevs;
    end
end



[s,in] = sort(allts(:,3));
allts   = allts(  in,:);
if nargout > 1
    allwave = allwave(in,:);
end