function yc = harmfit(x,y,nc,xc,err,isrobust)

% HARMFIT - Fit curve by harmonic functions
% yc = harmfit(x,y[,nc,xc,err]);
% x/y .. curve to be fitted
% yc  .. fitted curve
% nc  .. number of coefficients
% xc  .. abscissa to fit on (otherwise equal to x)

if nargin<6, isrobust=0; end
if (nargin<5)||isempty(err), err=ones(size(x)); end
if (nargin<4), 
    xc=x; 
else
    if length(xc)<2,
        isrobust=xc;
        xc=x;
    end
end
if (nargin<3)||(nc==0),
    nc=round(length(x)/30); % number of coefficients
end
xspan=max(x)-min(x);
xmi=min(x);
A=ones(length(x),nc*2+2)/2; %nc*(sin+cos)+offset+drift
A(:,2)=x(:)*3;
B=ones(length(xc),nc*2+2)/2;
B(:,2)=xc(:)*3;
for i=1:nc,
    A(:,i*2+1)=sin(2*i*pi*(x-xmi)/xspan);
    A(:,i*2+2)=cos(2*i*pi*(x-xmi)/xspan);
    B(:,i*2+1)=sin(2*i*pi*(xc(:)-xmi)/xspan);
    B(:,i*2+2)=cos(2*i*pi*(xc(:)-xmi)/xspan);
end
for i=0:isrobust*2,
    w=1./err(:);w(~isfinite(w))=0;
    WD=spdiags(w,0,length(w),length(w));
    coeff=(WD*A)\(y.*w);
    ww=irls(A*coeff-y);
    err(:)=err(:)./ww;
end
yc=B*coeff;
if nargout<1,
    plot(x,y,'b',xc,yc,'r-');
    rms(y,A*coeff)
end