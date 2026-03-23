function [rho,thk,R,S]=invfile1d(filename,nl,lam,thk,rho)

% INVFILE1D - 1D Inversion of File containing AB/2, MN/2 and rho_a
% [rho,thk,R] = invfile1d(filename,number_of_layers[,lambda]) OR
% [rho,thk,R] = invfile1d(matrix,number_of_layers[,lambda])
%               with matrix containing AB/2,MN/2 and rho_a columns

if ischar(filename), A=load(filename); else A=filename; end
ab2=A(:,1);mn2=A(:,2);R0=A(:,3);
if nargin<2, nl=4; end
if nargin<3, lam=0.5; end
aa=repmat(1:nl,2,1);inr=aa(:);inz=inr-1;inz(1:2)=[];
k=pi./(1./(ab2-mn2)-1./(ab2+mn2));
% err=100e-6*k/100/0.01;err*100
if nargin<5, rho=ones(1,nl)*median(R0); end
if nargin<4,
    thk=logspace(log10(min(ab2)),log10(max(ab2)),nl+1)/2;
    thk([1 end])=[];
end
maxz=max(ab2);minz=min(ab2)/3;
R=fwdschlum(rho,thk,ab2,mn2);
z=cumsum(thk);
if nargout<2,
    subplot(1,2,1);loglog(rho(inr),[minz z(inz) maxz],'b-');axis ij
    subplot(1,2,2);loglog(R0,ab2,'rx-',R,ab2,'bx-');axis ij;
end
it=0;rho1=rho;thk1=thk;tau=1;
rmserr=rms(R0,R);
fprintf('It %d: rms=%.1f%%\n',it,rmserr);
for it=1:20,
    S=senssch1d(rho,thk,ab2,mn2);
    dr=log(R0(:))-log(R(:));
    dm=(S'*S+lam*eye(size(S,2)))\(S'*dr);    
    rho(:)=rho(:).*exp(dm(1:nl));
    thk(:)=thk(:).*exp(dm(nl+1:end));
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
    end
    lam=lam*0.8;
end
if nargout>1,
    fprintf('rho=');fprintf('%d ',round(rho));
    fprintf('  dep=');fprintf('%.2f ',cumsum(thk));fprintf('\n');
end
if nargout<2, title(sprintf('RMS=%.1f%%',rms(R0,R))); end