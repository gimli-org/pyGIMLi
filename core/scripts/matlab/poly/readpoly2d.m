function Poly = readpoly2d(filename)

% READPOLY2D - Read triangle poly file to Poly struct
% Poly = readpoly2d(filename)

if nargin<1, filename='mesh.poly'; end
Poly=[];
fid=fopen(filename);
zeile=fgetl(fid);
aa=str2num(zeile);
Poly.dim=aa(2);
Poly.nnodes=aa(1);
ss='%d';for i=1:Poly.dim, ss=[ss '%f']; end
ss=[ss '%d'];
A=mytextscan(fid,ss,Poly.nnodes);
Poly.node=A{2};for i=2:Poly.dim+1, Poly.node=[Poly.node A{i+1}]; end
% Poly.nodemarker=A{Poly.dim+2};
zeile='';while isempty(zeile), zeile=fgetl(fid); end
aa=str2num(zeile);
Poly.nedges=aa(1);
ss='%d';for i=1:Poly.dim, ss=[ss '%f']; end
ss=[ss '%d'];
A=mytextscan(fid,ss,Poly.nedges);
Poly.edge=A{2}+1;for i=2:Poly.dim, Poly.edge=[Poly.edge A{i+1}+1]; end
Poly.edge=[Poly.edge A{Poly.dim+2}];
% Poly.egemarker=A{Poly.dim+2};
Poly.nholes=str2num(fgetl(fid));
Poly.hole=[];
for i=1:Poly.nholes,
   rr=str2num(fgetl(fid));
   Poly.hole(i,1:length(rr)-1)=rr(2:end);
end
Poly.nregions=str2num(fgetl(fid));
Poly.region=[];
for i=1:Poly.nregions,
   rr=str2num(fgetl(fid));
   Poly.region(i,1:length(rr)-1)=rr(2:end);
end
fclose(fid);
%%