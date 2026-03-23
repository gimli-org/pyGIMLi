function Poly = readpoly3d(filename)

% READPOLY3D - Read tetgen poly file to Poly struct
% Poly = readpoly3d(filename)

if nargin<1, filename='mesh.poly'; end
Poly=[];
fid=fopen(filename);
zeile=fgetl(fid);
aa=str2num(zeile);
Poly.nnodes=aa(1);
Poly.dim=aa(2);
ss='%d';for i=1:Poly.dim, ss=[ss '%f']; end
ss=[ss '%d'];
A=mytextscan(fid,ss,Poly.nnodes);
Poly.node=[A{2} A{3} A{4}];
Poly.nodemarker=A{5};
zeile='';while isempty(zeile), zeile=fgetl(fid); end
aa=str2num(zeile);
Poly.nfaces=aa(1);
Poly.facemarker=zeros(Poly.nfaces,1);
for i=1:Poly.nfaces,
   zeile='';while isempty(zeile), zeile=fgetl(fid); end
   aa=str2num(zeile);
   Poly.facemarker(i)=aa(3);
%    B=textscan(fid,'%d%f%f',aa(1));
%    Poly.face{i}=[B{2} B{3}];
   face={};
   for j=1:aa(1),
       aa=str2num(fgetl(fid));
       face{j}=aa(2:end)+1;
   end
   Poly.face{i}=face{1};%{:}
end
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
if nargout<1,
    for i=1:length(Poly.facemarker),
        if Poly.facemarker(i)<-2, % real face
            for j=1:length(Poly.face{i}),
                aa=Poly.face{i}{j};
                patch(Poly.node(aa,1),Poly.node(aa,2),Poly.node(aa,3),'b');
            end
        end
    end
    axis equal tight
end