function [res,thk,CHI2,resp]=invmaxmindata(data,nlay,coilsep,freq,lam,resL)

% [res,thk]=invmaxmindata(data,nlay,coilsep,freq,lam,resL)
% [res,thk]=invmaxmindata(data,(res,thk),coilsep,freq,lam,resL)

if nargin<6, resL=0; end
if nargin<5, lam=100; end
if nargin<4, freq=[110,220,440,880,1760,3520,7040,14080,28160,56320]; end
if nargin<3, coilsep=4; end
if nargin<2, nlay=4; end
zz=-1;
m0=0;nl=0;
if length(nlay)>1, % model given
    nl=fix((length(nlay)-1)/2)+1;
    m0=log(nlay);
    res=nlay(1:nl);
    thk=nlay(nl+1:nl*2-1);
    nlay=nl;
else
    res=ones(nlay,1)*10;
    thk=(1:nlay-1)'*3;
end
resp=maxminfwd(res,thk,coilsep,freq);
err=ones(length(freq)*2,1)*3;
err(10)=err(10)*10;
% err(18)=err(18)*10;
D=diag(1./err);
if nargout<1,
    showmaxmindata(data,resp,freq,3)
    subplot(131);draw1dmodel(res,thk);
end
for i=1:20,
    if nargout<1, pause(0.1); end
    S=maxminsens(res,thk,coilsep,freq,resL);
    DS=D*S;
    rhs=((D*(data(:)-resp(:)))'*DS)';
    if nl>0, rhs=rhs-([log(res-resL);log(thk)]-m0)*lam*0; end
    dm=(DS'*DS+eye(nlay*2-1)*lam)\rhs;
    if 1, %linesearch
        res1=(res-resL).*exp(dm(1:nlay))+resL;
        thk1=thk.*exp(dm(nlay+1:end));
        resp1=maxminfwd(res1,thk1,coilsep,freq);        
        taug=0.3;res03=(res-resL).*exp(dm(1:nlay)*taug)+resL;
        thk03=thk.*exp(dm(nlay+1:end)*taug);
        resp03=maxminfwd(res03,thk03,coilsep,freq);        
        taus=[0 taug 1];G=[1 0 0;1 taug taug^2;1 1 1];
        chis=[norm(D*(resp(:)-data(:)));norm(D*(resp03(:)-data(:)));norm(D*(resp1(:)-data(:)))];
        xx=G\chis;tau=min(-xx(2)/2/xx(3),1);
        if tau<=0.02, tau=0.1; end
    else
        tau=1;
    end
    res=(res-resL).*exp(dm(1:nlay)*tau)+resL;
    thk=thk.*exp(dm(nlay+1:end)*tau);
    resp=maxminfwd(res,thk,coilsep,freq);
    if nargout<1,
        showmaxmindata(data,resp,freq,3)
        subplot(131);draw1dmodel(res,thk);
    end
    lam=lam*0.8;
end
CHI2=chi2(data(:),resp(:),err,0);