function c=rndig(a,n)

% RNDIG - Round to n counting digits
% rndig(matrix,n)

if nargin<2, n=2; end
b=ones(size(a));
zn1=10^n;zn=zn1/10;
fi1=find((abs(a)<zn)&(a~=0)&isfinite(a));
fi2=find((abs(a)>zn1)&(a~=0)&isfinite(a));
while any(fi1)|any(fi2),
    b(fi1)=b(fi1)*10;
    a(fi1)=a(fi1)*10;
    b(fi2)=b(fi2)/10;
    a(fi2)=a(fi2)/10;
    fi1=find((abs(a)<zn)&(a~=0)&isfinite(a));
    fi2=find((abs(a)>zn1)&(a~=0)&isfinite(a));
end
c=round(a)./b;