function [A,coilsep,freq,x,y]=readmaxmindata(datafile)

% READMAXMINDATA - Read in MaxMin FDEM data from file
% [A,coilsep,freq,x,y] = readmaxmindata(datafile)

fid=fopen(datafile);
zeile=fgetl(fid);
coilsep=sscanf(strrep(zeile,':',' '),'%*s%*s%f%*s');
zeile=fgetl(fid);
fstr=strrep(zeile(strfind(zeile,':')+1:strfind(zeile,'Hz')-1),',',' ');
freq=str2num(fstr);
zstar=sqrt(30./freq)*500;
fclose(fid);
A=textread(datafile,'','headerlines',4);
y=A(:,1);
x=A(:,2);
A(:,1:2)=[];
