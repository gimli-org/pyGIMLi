function N=readunifile(fname,dim)

% READUNIFILE - Read data file in unified data format
% N = readunifile(filename)
% Format: (# character can be used to comment)
% number_of_Electrodes
% # x y h d # token string for meaning of columns
% x_el1 (y_el1) z_el1
% ...
% x_eln (y_eln) z_eln
% number_of_datapoints
% # a b m n u i # possible token string for meaning of columns
% A_1 B_1 M_1 N_1 R_1 (Err_1)
% ...
% A_n B_n M_n N_n R_n (Err_n)
% (= Electrode numbers, 0=infinity)

if nargin<2, dim=0; end
N.elec=[];
fid=fopen(fname,'r');
if fid<0, error('File not found!'); end
zeile=destrip(fgetl(fid));
ne=sscanf(zeile,'%d\n',1);
zeile=fgetl(fid);formstr='';
if (zeile(1)=='#')&&(strfind(lower(zeile),'x')||strfind(lower(zeile),'y')),
    ix=0;iy=0;iz=0;ih=0;id=0;xmul=1;ymul=1;zmul=1;hmul=1;dmul=1;
    zeile(1)='';zeile=destrip(zeile);i=0;formstr='';
    while ~isempty(zeile),
        i=i+1;utok='';
        [tok,zeile]=strtok(zeile);
        if isempty(tok)||isequal(tok,'position'), break; end
        fis=strfind(tok,'/');
        if ~isempty(fis), % physical unit found
            fis=fis(1);utok=tok(fis+1:end);tok=tok(1:fis-1);
        end
        mul=1;
        switch utok,
            case 'm', mul=1;
            case 'cm', mul=0.01;
            case 'in', mul=0.0254;
            case 'mm', mul=1e-3;
            case 'km', mul=1e3;
            otherwise, if ~isempty(utok), fprintf('Unit %s not recognized\n',utok); end
        end
        fprintf('%s ',tok);
        switch lower(tok),
            case 'x', ix=i;xmul=mul;formstr=[formstr '%f'];
            case 'y', iy=i;ymul=mul;formstr=[formstr '%f'];
            case 'z', iz=i;zmul=mul;formstr=[formstr '%f'];
            case 'h', ih=i;hmul=mul;formstr=[formstr '%f'];
            case 'd', id=i;dmul=mul;formstr=[formstr '%f'];
            otherwise, fprintf('Token %s not recognized\n',tok);
        end
    end
    A=mytextscan(fid,formstr,ne,'commentstyle','#');
    % convention that h/d have to be used instead of z for future storage
    % if (iz>0)&&(ih==0)&&(id==0), id=iz;ih=iz;dmul=zmul;hmul=zmul;iz=0; end
    if ix, N.x=A{ix}; end
    if iy, N.y=A{iy}; end
    if iz, N.z=A{iz}; end
    if ih, N.h=A{ih}; end
    if id, N.d=A{id}; end
    if ix, N.elec=[N.elec N.x]; end
    if iy, N.elec=[N.elec N.y]; end
    if iz, N.elec=[N.elec N.z]; end
    if dim==0, dim=max((ix>0)+(iy>0)+(iz>0),2); end
    if (dim==2)&(size(N.elec,2)>1), % x and y present
        if length(unique(N.elec(:,1)))==1, N.elec(:,1)=[];
        elseif length(unique(N.elec(:,2)))==1, N.elec(:,2)=[];
        else N.elec=(0:length(N.x)-1)';
        end
    end
    %         if (ih>0)&&(iz>0), N.elec=[N.elec N.h-abs(N.d)];
    %         elseif (iz>0), N.elec=[N.elec N.z];
    %         elseif (id>0), N.elec=[N.elec abs(N.d)];
    %         else N.elec(:,dim)=0;
    %         end
    if (ih>0)&&(iz>0), N.elec(:,dim)=N.h-abs(N.d);
    elseif (iz>0), N.elec(:,dim)=N.z;
    elseif (id>0), N.elec(:,dim)=abs(N.d);
    else N.elec(:,dim)=0;
    end
    if size(N.elec,2)>dim, N.elec(:,dim+1:end)=[]; end
else
    for n=1:ne,
        if n>1, zeile=destrip(fgetl(fid)); end
        while isempty(zeile), zeile=destrip(fgetl(fid)); end
        el=str2num(zeile);
        N.elec(n,1:length(el))=el;
    end
    if size(N.elec,2)<2, N.elec(:,2)=0; end
    if size(N.elec,2)<3, N.elec(:,3)=0; end
end
zeile='';
while isempty(zeile), zeile=destrip(fgetl(fid)); end
nm=sscanf(zeile,'%d\n',1);
zeile='';
while isempty(zeile), zeile=fgetl(fid); end
lz=lower(zeile);
if (zeile(1)=='#')&&(any(strfind(lz,'a'))||...
        any(strfind(lz,'b'))||any(strfind(lz,'s'))),
    ia=0;ib=0;im=0;in=0;ir=0; % A B M N Rhoa
    ii=0;iu=0;ik=0; % current voltage geometric factor
    it=0;ig=0;is=0; % time geophone shot
    ierr=0;iip=0;ifr=0;ieip=0; %err ip frequency iperror
    emul=1;imul=1;umul=1;tmul=1;irho=0;isp=0;
    zeile(1)='';zeile=destrip(zeile);i=0;formstr='';
    while ~isempty(zeile),
        i=i+1;utok='';
        [tok,zeile]=strtok(zeile);
        fis=strfind(tok,'/');
        if ~isempty(fis), % physical unit found
            fis=fis(1);utok=tok(fis+1:end);tok=tok(1:fis-1);
        end
        if ~isempty(tok),
            fprintf('%s ',tok)
            switch lower(tok),
                case {'a','c1'}, ia=i;formstr=[formstr '%d'];
                case {'b','c2'}, ib=i;formstr=[formstr '%d'];
                case {'m','p1'}, im=i;formstr=[formstr '%d'];
                case {'n','p2'}, in=i;formstr=[formstr '%d'];
                case {'g','geophone'}, ig=i;formstr=[formstr '%d'];
                case {'s','shot'}, is=i;formstr=[formstr '%d'];
                case {'rhoa','ra','rho_a'}, ir=i;formstr=[formstr '%f'];
                case {'rho','r'}, irho=i;formstr=[formstr '%f'];
                case {'err','error','std'},
                    ierr=i;formstr=[formstr '%f'];
                    if isequal(utok,'%'), emul=0.01; end
                case 'ip', iip=i;formstr=[formstr '%f'];
                case 'iperr', ieip=i;formstr=[formstr '%f'];
                case 'sp', isp=i;formstr=[formstr '%f'];
                case 'f', ifr=i;formstr=[formstr '%f'];
                case {'i','cur','current'}, ii=i;formstr=[formstr '%f'];
                    if isequal(utok,'mA'), imul=1e-3; end
                    if isequal(utok,'uA'), imul=1e-6; end
                    if isequal(utok,'nA'), imul=1e-9; end
                    if isequal(utok,'kA'), imul=1e+3; end
                case {'u','v','volt','voltage'}, iu=i;formstr=[formstr '%f'];
                    if isequal(utok,'mV'), umul=1e-3; end
                    if isequal(utok,'uV'), umul=1e-6; end
                    if isequal(utok,'nV'), umul=1e-9; end
                    if isequal(utok,'kV'), umul=1e+3; end
                case {'k','g'}, ik=i;formstr=[formstr '%f'];
                case {'t','topo','tt','traveltime'}, it=i;formstr=[formstr '%f'];
                    if isequal(utok,'ms'), tmul=1e-3; end
                    if isequal(utok,'us'), tmul=1e-6; end
                    if isequal(utok,'ns'), tmul=1e-9; end
                otherwise, %unknown token
                    formstr=[formstr '%*s'];
                    i=i-1;
                    fprintf('(%s not found) ',tok);
            end
        end
    end
    fprintf('found\n');
    %     [ia ib im in ir ierr iip ii iu],formstr
    A=mytextscan(fid,formstr,nm,'commentstyle','#');
    if ia, N.a=A{ia}; end
    if ib, N.b=A{ib}; end
    if im, N.m=A{im}; end
    if in, N.n=A{in}; end
    if ig, N.g=A{ig}; end
    if is, N.s=A{is}; end
    if ifr, N.f=A{ifr}; end
    if ir, N.r=A{ir}; end
    if irho, N.rho=A{irho}; end
    if iip, N.ip=A{iip}; end
    if ieip, N.iperr=A{ieip}; end
    if isp, N.sp=A{isp}; end
    if ierr, N.err=A{ierr}*emul; end
    if ii, N.i=A{ii}*imul; end
    if iu, N.u=A{iu}*umul; end
    if ik, N.k=A{ik}; end
    if it, N.t=A{it}*tmul; end
    if isfield(N,'a')|isfield(N,'b')|isfield(N,'m')|isfield(N,'n'), % DC files
        if ~isfield(N,'b'), N.b=zeros(size(N.a)); end
        if ~isfield(N,'n'), N.n=zeros(size(N.m)); end
    end
else
    N.a=zeros(nm,1);N.b=N.a;N.m=N.a;N.n=N.a;N.r=N.a;N.err=N.a;N.ip=N.a;
    zeile=destrip(zeile);
    for n=1:nm,
        if n>1, zeile=destrip(fgetl(fid)); end
        while isempty(zeile), zeile=destrip(fgetl(fid)); end
        mess=str2num(zeile);
        if length(mess)<5, break; end
        N.a(n)=mess(1);N.b(n)=mess(2);N.m(n)=mess(3);N.n(n)=mess(4);
        if length(mess)>4, N.r(n)=mess(5); end
        if length(mess)>5, N.err(n)=mess(6); end
        if length(mess)>6, N.ip(n)=mess(7); end
    end
    if max(N.ip)<=0, N=rmfield(N,'ip'); end
end
zeile='';
while isempty(zeile)&&(~feof(fid)), zeile=destrip(fgetl(fid)); end
if ~feof(fid),
    ntopo=str2num(zeile);
    N.topo=zeros(ntopo,2);
    for i=1:ntopo,
        zeile='';
        while isempty(zeile), zeile=destrip(fgetl(fid)); end
        xz=str2num(zeile);
        if length(xz)>1, N.topo(i,:)=xz(1:2); end
    end
end
fclose(fid);
if isfield(N,'a'),
    nm=min([length(N.a) length(N.b) length(N.m) length(N.n)]);
end
fn=fieldnames(N);
for i=1:length(fn),
    fie=getfield(N,fn{i});
    %             if (min(size(fie))==1)&&(length(fie)>nm), fie(nm+1:end)=[];N=setfield(N,fn{i},fie); end
    if ~ismember(fn{i},{'x','y','z'})&&(min(size(fie))==1)&&(length(fie)>nm),
        fie(nm+1:end)=[];N=setfield(N,fn{i},fie); end
end
% if isfield(N,'err')&&(min(N.err)<=0),
%     messg('Found nonpositive errors. Discarding error values.');
%     N=rmfield(N,'err');
% end

% if ~isfield(N,'k')||(length(N.k)<length(N.a)), N.k=getkonf(N); end

% if ~isfield(N,'r'), % no apparent resistivities
%     if ~isfield(N,'rho')&&isfield(N,'u')&&isfield(N,'i'), N.rho=N.u./N.i; end
%     if isfield(N,'rho'), N.r=N.rho.*N.k; end
% end
display(sprintf('%s: %d Measurements with %d Electrodes',fname,nm,ne));

function zeile=destrip(zeile)
% strip string from comments (with # character)
aa=strfind(zeile,'#');
if ~isempty(aa), zeile=zeile(1:aa(1)-1); end
