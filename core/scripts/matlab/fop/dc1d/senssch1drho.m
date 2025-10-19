function S = senssch1drho(rho,thk,ab2,mn2,islin)

% SENSSCH1D - Sensitivity for schlumberger 1d
% S = senssch1d(rho,thk,ab2,mn2)

if nargin<5, islin=0; end
if ~isequal(size(ab2),size(mn2)), error('AB/2 must equal MN/2'); end
if length(rho)~=length(thk)+1, error('Resistivity must equal thicknesses+1'); end
S=zeros(length(ab2),length(rho));
R0=fwdschlum(rho,thk,ab2,mn2);
fak=1.05;
for i=1:length(rho),
    rho1=rho;rho1(i)=rho(i)*fak;
    R=fwdschlum(rho1,thk,ab2,mn2);
    if islin,
        S(:,i)=(R(:)-R0(:))/((fak-1)*rho(i));
    else
        S(:,i)=(log(R(:))-log(R0(:)))/log(fak);        
    end
end