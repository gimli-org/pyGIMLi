function errorbarxy(x,y,lx,ly,ux,uy,linecol,errorcol)

% ERRORBARXY - Error bars for both x and y
% errorbarxy(x,y,lx,ly,ux,uy)

%This function allows the user to plot the graph of x against y, along with both x and y errorbars.
%For the x and y errors it is possible to input both lower (lx and ly)  and upper  (ux and uy) values for the
%errors at a particular point.  If the upper values are not specified then the program assumes the errors 
%are symmetrical and use the lower values.  it is also possible to specify the plot line colour, marker, and 
%linestyle using the standard 'plot' command notation in the input variable 'linecol'.  Also the line colour for 
%the errobars can be specified in the variable 'errorcol'.  It is important to note that if these colour options 
%are to be used and any of the error limit vectors are empty then they should not be excluded, but presented 
%in a [] form signifying an empty vector.
%
%James Rooney,  17 October 2003


if exist('linecol','var')==0 | isempty(linecol)
    linecol='bx-';
end


if exist('errorcol','var')==0 | isempty(errorcol)
    errorcol='r';
end

if nargin<4, ly=zeros(size(y)); end

plot(x,y,linecol)
hold on

xw=ones(size(x))*(max(x)-min(x))/100;
yw=ones(size(x))*(max(y)-min(y))/100;
yw=yw/10;

lye=exist('ly','var');
lxe=exist('lx','var');
uye=exist('uy','var');
uxe=exist('ux','var');

if lye+lxe+uye+uxe==0 | isempty(lx) & isempty(ux) & isempty(ly) & isempty(uy)
    return
end

if uye==0 | isempty(uy)
    uy=ly;
end

if uxe==0 | isempty(ux)
    ux=lx;
end

for t=1:length(x)
    
if ~isempty(ux)
    %x errorbars
        line([x(t)-lx(t) x(t)+ux(t)],[y(t) y(t)],'color',errorcol)
    line([x(t)-lx(t) x(t)-lx(t)],[y(t)-yw(t) y(t)+yw(t)],'color',errorcol)    
        line([x(t)+ux(t) x(t)+ux(t)],[y(t)-yw(t) y(t)+yw(t)],'color',errorcol)    
    end
    
    if ~isempty(uy)
    %y errorbars
    line([x(t) x(t)],[y(t)-ly(t) y(t)+uy(t)],'color',errorcol)
    line([x(t)-xw(t) x(t)+xw(t)],[y(t)-ly(t) y(t)-ly(t)],'color',errorcol)    
        line([x(t)-xw(t) x(t)+xw(t)],[y(t)+uy(t) y(t)+uy(t)],'color',errorcol) 
    end    
    end
    hold off
    