function [ temp, date, time ] = GetData( file )
fileID = fopen(file,'r');
formatSpec = '%c' ;
A  = fscanf(fileID,formatSpec);
B = strsplit(A,';');

xtemp = 3:3:length(B);
temp = B(xtemp);

xdate = 1:3:length(B);
date = B(xdate);

xtime = 2:3:length(B);
time = B(xtime);

end

