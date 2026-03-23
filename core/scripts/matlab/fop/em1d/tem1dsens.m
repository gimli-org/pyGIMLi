function S = tem1dsens(res,thk,T,txarea,rxarea)

nlay=min(length(res),length(thk)+1);
S=zeros(length(T),nlay*2-1);
rhoa=tem1drhoa(res,thk,T,txarea,rxarea);
fak=1.05;
for i=1:nlay,
    res1=res;
    res1(i)=res(i)*fak;
    rhoa1=tem1drhoa(res1,thk,T,txarea,rxarea);
    S(:,i)=(log(rhoa1)-log(rhoa))/log(fak);
end
for i=1:nlay-1,
    thk1=thk;
    thk1(i)=thk(i)*fak;
    rhoa1=tem1drhoa(res,thk1,T,txarea,rxarea);
    S(:,nlay+i)=(log(rhoa1)-log(rhoa))/log(fak);
end
