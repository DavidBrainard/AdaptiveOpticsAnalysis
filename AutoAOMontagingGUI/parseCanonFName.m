function [ eyeSide, LocXY ] = parseCanonFName( filenameStr )



%% File structure:
% ???-YYYMMDD???????-SUBID-??-SEX-EYE--XXXYYY-SCANTYPE-SCANNUM-EXPORTINFO
% 206-20140730132907-11002-00-F-OS--063043-P0-0002-Preprocessed20160630104614

tokenizedStr = textscan(filenameStr,'%s','Delimiter','-');
tokenizedStr = tokenizedStr{1};
%% FOR CANON FILENAMES POST Nov 13, 2012
eyeSide = tokenizedStr{6};

locStr = tokenizedStr{8};

% First part is the x location
LocXY(1) = 0.5 * (str2double(locStr(1:2)) - 6);
% Second part is the y location
LocXY(2) = 0.5 * (str2double(locStr(4:5)) - 6);


end