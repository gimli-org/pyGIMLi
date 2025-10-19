function [rhoa,T,txarea,rxarea,V]=readtemfile(temfile)

fid=fopen(temfile);
isdata=0;
A=[];
while 1,
   zeile=fgetl(fid);
   if isempty(zeile), isdata=0; end
   if strfind(zeile,'<END FILE>'), break; end
   if isdata,
       A=[A;str2num(zeile)];
   end
   if strfind(zeile,'TXAREA'), txarea=sscanf(zeile,'%*s%f%*s'); end
   if strfind(zeile,'RXAREA'), rxarea=sscanf(zeile,'%*s%f%*s'); end
   if strfind(zeile,'# Time'), isdata=1; end
end
fclose(fid);

T=A(:,1)/1000;
V=A(:,2)/1e9;
%%
mu0=4e-7*pi;
a=sqrt(txarea/pi);
fak=a^(4/3)*rxarea^(2/3)*mu0^(5/3)/20^(2/3)/pi^(1/3);
rhoa=fak./(T.^(5/3).*V.^(2/3));
