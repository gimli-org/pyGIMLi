function [x,y,z,Zone]=gps2xyz(lat,lon,h,lcm)

% GPS2XYZ - converts gps output (5412.09) to xyz
%           based on the geodetic toolbox by Michael R. Craymer 
% [x,y,z] = gps2xyz(lat,lon[,h])

if nargin<3, h=zeros(size(lat)); end
a=6378137.0;finv=298.257223563; %WGS84
e2=1-(1-1/finv)^2; 

de=floor(lat/100);
mi=floor(lat-de*100);
se=(lat-mi-de*100)*60;
nlat=dms2rad([de mi se]);
de=floor(lon/100);
mi=floor(lon-de*100);
se=(lon-mi-de*100)*60;
nlon=dms2rad([de mi se]);
% [x,y,z]=ell2xyz(nlat,nlon,h,a,e2);
if nargin>3,
    [x,y,z,Zone]=ell2utm(nlat,nlon,a,e2,lcm);
else
    [x,y,z,Zone]=ell2utm(nlat,nlon,a,e2);
end
% x=-x;%y=-y;

function [N,E,h,Zone]=ell2utm(lat,lon,a,e2,lcm)
% ELL2UTM  Converts ellipsoidal coordinates to UTM.
%   UTM northing and easting coordinates in a 6 degree
%   system.  Zones begin with zone 1 at longitude 180E
%   to 186E and increase eastward.  Formulae from E.J.
%   Krakiwsky, "Conformal Map Projections in Geodesy",
%   Dept. Surveying Engineering Lecture Notes No. 37,
%   University of New Brunswick, Fredericton, N.B.
%   Vectorized.
% Version: 31 Mar 2005
% Useage:  [N,E,Zone]=ell2utm(lat,lon,a,e2,lcm)
% Input:   lat - vector of latitudes (rad)
%          lon - vector of longitudes (rad)
%          a   - major semi-axis of ref. ell (m)
%          e2  - eccentricity squared of ref. ell.
%          lcm - optional central meridian (default =
%                standard UTM den'n)
% Output:  N   - vector of UTM northings (m)
%          E   - vector of UTM eastings (m)
%          Zone- vector of UTM zones

ko=0.9996;      % Scale factor
if lat>=0
  No=0;         % False northing (north)
else
  No=10000000;  % False northing (south)
end
Eo=500000;      % False easting

if nargin==5
  Zone=zeros(size(lat));
else
  Zone=floor((rad2deg(lon)-180)/6)+1;
  Zone=Zone+(Zone<0)*60-(Zone>60)*60;
  lcm=deg2rad(Zone*6-183);
end

lam=lon-lcm;
lam=lam-(lam>=pi)*(2*pi);
  
%fprintf('\nZones\n');
%fprintf('%3d\n',Zone');
%fprintf('\nCentral Meridians\n');
%fprintf('%3d %2d %9.6f\n',rad2dms(lcm)');
%fprintf('\nLongitudes wrt Central Meridian\n');
%fprintf('%3d %2d %9.6f\n',rad2dms(lam)');

f=1-sqrt(1-e2);
RN=a./(1-e2*sin(lat).^2).^0.5;
RM=a*(1-e2)./(1-e2*sin(lat).^2).^1.5;
t=tan(lat);
h=sqrt(e2*cos(lat).^2/(1-e2));
n=f/(2-f);

a0=1+n^2/4+n^4/64;
a2=1.5*(n-n^3/8);
a4=15/16*(n^2-n^4/4);
a6=35/48*n^3;
a8=315/512*n^4;

s=a/(1+n)*(a0*lat-a2*sin(2*lat)+a4*sin(4*lat)- ...
  a6*sin(6*lat)+a8*sin(8*lat));

E1=lam .* cos(lat);
E2=lam.^3 .* cos(lat).^3/6 .* (1-t.^2+h.^2);
E3=lam.^5 .* cos(lat).^5/120 .* ...
    (5-18*t.^2+t.^4+14*h.^2-58*t.^2 .*h.^2+13*h.^4+...
     4*h.^6-64*t.^2 .*h.^4-24*t.^2 .*h.^6);
E4=lam.^7 .*cos(lat).^7/5040 .* ...
    (61-479*t.^2+179*t.^4-t.^6);
E=Eo + ko*RN.*(E1 + E2 + E3 + E4);

N1=lam.^2/2 .* sin(lat) .* cos(lat);
N2=lam.^4/24 .* sin(lat) .* cos(lat).^3 .* ...
    (5-t.^2+9*h.^2+4*h.^4);
N3=lam.^6/720 .* sin(lat) .* cos(lat).^5 .* ...
    (61-58*t.^2+t.^4+270*h.^2-...
     330*t.^2 .*h.^2+445*h.^4+...
     324*h.^6-680*t.^2 .*h.^4+...
     88*h.^8-600*t.^2 .*h.^6-...
     192*t.^2 .*h.^8);
N4=lam.^8/40320 .* sin(lat) .* cos(lat).^7 .* ...
   (1385-311*t.^2+543*t.^4-t.^6);
N=No + ko*RN.*(s./RN + N1 + N2 + N3 + N4); 

function deg=rad2deg(rad)
% RAD2DEG  Converts radians to decimal degrees. Vectorized.
% Version: 8 Mar 00
% Useage:  deg=rad2deg(rad)
% Input:   rad - vector of angles in radians
% Output:  deg - vector of angles in decimal degrees
deg=rad.*180./pi;
%ind=(deg<0);
%deg(ind)=deg(ind)+360; 

function rad=deg2rad(deg)
% DEG2RAD  Converts decimal degrees to radians. Vectorized.
% Version: 18 Jan 96
% Useage:  rad=deg2rad(deg)
% Input:   deg - vector of angles in decimal degrees
% Output:  rad - vector of angles in radians
rad=deg.*pi./180; 

function [x,y,z]=ell2xyz(lat,lon,h,a,e2)
% ELL2XYZ  Converts ellipsoidal coordinates to cartesian.
%   Vectorized.
% Version: 18 Jan 96
% Useage:  [x,y,z]=ell2xyz(lat,lon,h,a,e2)
% Input:   lat - vector of ellipsoidal latitudes (radians)
%          lon - vector of ellipsoidal E longitudes (radians)
%          h   - vector of ellipsoidal heights (m)
%          a   - ref. ellipsoid major semi-axis (m)
%          e2  - ref. ellipsoid eccentricity squared
% Output:  x \
%          y  > vectors of cartesian coordinates in CT system (m)
%          z /
v=a./sqrt(1-e2*sin(lat).*sin(lat));
x=(v+h).*cos(lat).*cos(lon);
y=(v+h).*cos(lat).*sin(lon);
z=(v.*(1-e2)+h).*sin(lat); 

function rad=dms2rad(dms)
% DMS2RAD  Converts degrees-minutes-seconds to radians.
%   Vectorized.
% Version: 12 Mar 00
% Useage:  rad=dms2rad(dms)
% Input:   dms - [d m s] array of angles in deg-min-sec, where
%                d = vector of degrees
%                m = vector of minutes
%                s = vector of seconds
% Output: rad - vector of angles in radians
d=dms(:,1);
m=dms(:,2);
s=dms(:,3);
dec=abs(d)+abs(m)./60+abs(s)./3600;
rad=dec.*pi./180;
ind=(d<0 | m<0 | s<0);
rad(ind)=-rad(ind); 