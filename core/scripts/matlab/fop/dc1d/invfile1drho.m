function [rho,R]=invfile1drho(filename,thk,lam,rho)

% INVFILE1DRHO - 1D Inversion of File containing AB/2, MN/2 and rho_a
% [rho,R] = invfile1drho(filename,thicknesses[,lambda]) OR
% [rho,R] = invfile1drho(matrix,thicknesses[,lambda])
%               with matrix containing AB/2,MN/2 and rho_a columns

if ischar(filename), A=load(filename); else A=filename; end
ab2=A(:,1);mn2=A(:,2);R0=A(:,3);
if nargin<2, nl=4; end
if nargin<3, lam=0.5; end
nl=length(thk)+1;
if nargin<4, rho=ones(1,nl)*median(R0); end
aa=repmat(1:nl,2,1);inr=aa(:);inz=inr-1;inz(1:2)=[];
k=pi./(1./(ab2-mn2)-1./(ab2+mn2));
% err=100e-6*k/100/0.01;err*100
maxz=max(ab2);minz=min(ab2)/3;
R=fwdschlum(rho,thk,ab2,mn2);
z=cumsum(thk);
if nargout<2,
    subplot(1,2,1);loglog(rho(inr),[minz z(inz) maxz],'b-');axis ij
    subplot(1,2,2);loglog(R0,ab2,'rx-',R,ab2,'bx-');axis ij;
    legend('measured','calculated',2);
end
it=0;rho1=rho;tau=1;
rmserr=rms(R0,R);
fprintf('It %d: rms=%.1f%%\n',it,rmserr);
len=length(rho);
if 1,
    C=eye(len);
else
    C=zeros(len-1,len);
    for i=1:len-1, C(i,i:i+1)=[-1 1]; end
end
CTC=C'*C;
for it=1:10,
    S=senssch1drho(rho,thk,ab2,mn2);
    dr=log(R0(:))-log(R(:));
    dm=(S'*S+lam*CTC)\(S'*dr-lam*CTC*log(rho(:)));
    for t=1:10,
        rho1=rho(:).*exp(dm*t/10);
        R1=fwdschlum(rho1,thk,ab2,mn2);
        rmst(t)=rms(R0,R1);
    end
    [mi,itau]=min(rmst);tau=itau/10;
    rho(:)=rho(:).*exp(dm*tau);
    R=fwdschlum(rho,thk,ab2,mn2);
    oldrms=rmserr;rmserr=rms(R0,R);
    fprintf('It %d: rms=%.1f%%\n',it,rms(R0,R));
    if rmserr/oldrms>0.9, break; end
    z=cumsum(thk);
    if nargout<2, 
        subplot(1,2,1);loglog(rho(inr),[minz z(inz) maxz],'b-');axis ij
        xlabel('rho in Ohmm');ylabel('z in m');grid on;
        subplot(1,2,2);loglog(R0,ab2,'rx-',R,ab2,'bx-');axis ij;
        xlabel('rhoa in Ohmm');ylabel('AB/2 in m');grid on;
        legend('measured','calculated');
    end
end
if nargout>1,
    fprintf('rho=');fprintf('%d ',round(rho));
    fprintf('  dep=');fprintf('%.2f ',cumsum(thk));fprintf('\n');
end
if nargout<2, title(sprintf('RMS=%.1f%%',rms(R0,R))); end

function S = senssch1drho(rho,thk,ab2,mn2)

% SENSSCH1D - Sensitivity for schlumberger 1d
% S = senssch1d(rho,thk,ab2,mn2)

if ~isequal(size(ab2),size(mn2)), error('AB/2 must equal MN/2'); end
if length(rho)~=length(thk)+1, error('Resistivity must equal thicknesses+1'); end
S=zeros(length(ab2),length(rho));
R0=fwdschlum(rho,thk,ab2,mn2);
fak=1.05;
for i=1:length(rho),
    rho1=rho;rho1(i)=rho(i)*fak;
    R=fwdschlum(rho1,thk,ab2,mn2);
    S(:,i)=(log(R(:))-log(R0(:)))/log(fak);
end
	
