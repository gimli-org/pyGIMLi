function show2dmrsimage(data,meas,field)

% SHOW2DMRSIMAGE - Show 2d MRS data as subplot curves
% show2smrsimage(data,meas,field)

if nargin<3, 
    field=abs(data.V_mes);
end
cla;
Nqmax=max(diff(data.idx_q));
A=ones(Nqmax,data.N_cfg)*NaN;
for i=1:data.N_cfg,
   idx=data.idx_q(i):data.idx_q(i+1)-1;
   [qq,sq]=sort(data.q_val(idx));
   A(1:length(idx),i)=field(idx(sq));
end
AL=1-isnan(A);
imagesc(A);
% caxis([0 0.3]);
colorbar;
alpha(AL);
% colormap(flipud(jet));
for j=1:data.N_cfg, 
    xtl{j}=[num2str(meas.ntx(j)) '-' num2str(meas.nrx(j))]; 
end
set(gca,'XTick',1:size(A,1),'XTickLabel',xtl);
