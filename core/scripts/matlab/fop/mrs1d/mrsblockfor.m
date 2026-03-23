function [u,wcvec]=mrsblockfor(Kerneldata,wc,thk)

% MRSBLOCKFOR - MRS Block forward calculation
% mrsblockfor(Kerneldata,wc,thk)
% Kerneldata..struct of K (kernel) and z or dz

if isfield(Kerneldata.model,'z'), %
    z=[0;Kerneldata.model.z(:)];
elseif 1, % 
    z=[0;cumsum(Kerneldata.model.dz(1:end-1))'];
else % old version
    z=(0:size(Kerneldata.K,1)-1).*Kerneldata.model.dz;
end
zhk=cumsum(thk);
wcvec=ones(size(Kerneldata.K,1),1)*wc(end);
iz1=0;
for i=1:length(thk),
    iz2=find(z<zhk(i),1,'last');
    wcvec(iz1+1:iz2)=wc(i);
    if iz2>=length(z), break; end
    wcvec(iz2)=((zhk(i)-z(iz2))*wc(i)+(z(iz2+1)-zhk(i))*wc(i+1))/(z(iz2+1)-z(iz2));
    iz1=iz2;
end
u=abs(Kerneldata.K'*wcvec);