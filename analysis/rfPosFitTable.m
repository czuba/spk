function rftable = rfPosFitTable(rfStruct)

% extract fit:  x, y, ecct, fwhm, r2
tableDat = [rfStruct.fit(:,3), rfStruct.fit(:,4), hypot(rfStruct.fit(:,3), rfStruct.fit(:,4)), sqrt(rfStruct.fit(:,5))*2.355, rfStruct.r2];

rftable = array2table(tableDat, 'VariableNames',{'xDeg','yDeg','Ecct','FWHM','r2'},'RowNames',compose('%d',1:size(rfStruct.fit,1)));