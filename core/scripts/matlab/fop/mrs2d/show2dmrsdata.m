function show2dmrsdata(data,meas,field,field2)

% SHOW2DMRSDATA - Show 2d MRS data as subplot curves
% show2smrsdata(data,meas[,field[,field2]])

if (nargin<3)||isempty(field), 
    field=abs(data.V_mes);
end
Ntx=max(meas.ntx);
Nrx=max(meas.nrx);
qmi=min(data.q_val);
qmi=min(qmi,rndig(qmi,1)); % not too much zeros
qma=max(data.q_val);
fmi=min(field);
fma=max(field);
clf;
for i=1:data.N_cfg,
   subplot(Ntx,Nrx,data.mes_idx(i));
   idx=data.idx_q(i):data.idx_q(i+1)-1;
   semilogy(field(idx),data.q_val(idx),'bx-');
   if nargin>3,
       hold on;
       semilogy(field2(idx),data.q_val(idx),'r+-');
       hold off;
   end
   axis ij;grid on;
   xlim([fmi fma]);ylim([qmi qma]);
   yt=get(gca,'YTick');
   if length(yt)<3,
      yt=[rndig(qmi,1);yt(:);qma];      
      ytl=num2strcell(rndig(yt,2));
      set(gca,'YTick',yt,'YTickLabel',ytl);
   end
   
   text(fmi,qmi,['Tx' num2str(meas.ntx(i)) 'Rx' num2str(meas.nrx(i))],...
        'HorizontalAlignment','left','VerticalAlignment','bottom');
end