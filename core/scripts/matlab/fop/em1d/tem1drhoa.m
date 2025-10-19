function rhoa=tem1drhoa(res,thk,T,txarea,rxarea)

% TEM1DRHOA - Calculate TEM 1d apparent resistivity
% rhoa=tem1drhoa(res,d,T,txarea,rxarea)

if nargin<5, rxarea=txarea; end
mu0=4e-7*pi;
a=sqrt(txarea/pi);
fak=a^(4/3)*rxarea^(2/3)*mu0^(5/3)/20^(2/3)/pi^(1/3);
[HZP,HZ,HR,EPHI0,RR,T2]=vmd_t([1,1,1,1],min(T),max(T),res,thk,a,a,0,rxarea);
EPHI=exp(interp1(log(T2),log(EPHI0),log(T)));
rhoa=fak./(T.^(5/3).*(EPHI*2*pi*a).^(2/3));
