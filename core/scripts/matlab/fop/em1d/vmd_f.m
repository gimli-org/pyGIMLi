function [fields,rr]=vmd_f(iout,f,rho,d,ze,zs,rmin,nr,tm)
% [FIELDS,RR]=VMD_F(IOUT,F,RHO,D,ZE,ZS,RMIN,NR,TM)
%
% Berechnung der Felder eines zeitharmonischen Vertikalen
% Magnetischen Dipols (VMD) ueber einem geschichteten Halbraum
% mittels Hankeltransformation.
% Verwendet werden Zylinderkoordinaten mit z positiv nach
% unten.
% Der VMD befindet sich in r=0 und z=zs, zs<0. Die
% Horizontalkomponenten der Aufpunkte RR sind logarithmisch
% aequidistant verteilt und werden aus RMIN und NR berechnet. Die
% Felder werden fuer z=ze bestimmt.
% Die verwendeten Filterkoeffizienten sind fuer 10 Stuetzstellen
% pro Dekade vorgegeben.
% Der geschichtete Halbraum wird durch Angabe von NL
% spezifischen Widerstaenden und NL-1 Maechtigkeiten beschrieben.
%
% Input:
% IOUT:  Auswahl der zu berechnenden Felder
%        IOUT(3)
%        Vorschrift: [HZ,HR,EPHI].*IOUT
%        Beispiele: IOUT=[1,0,0]: Berechnung von HZ
%                   IOUT=[0,1,1]: Berechnung von HR und EPHI
% F:     Frequenz in Hertz
% RHO:   Spezifische Widerstaende in Ohm*m, RHO(NL)
% D:     Maechtigkeiten in m, D(NL-1), [] fuer hom. Halbraum
% ZE:    z-Koordinate des Empfaengers in Meter (default: ZE=0)
% ZS:    z-Koordinate des Senders in Meter (default: ZS=0)
% RMIN:  Kleinster Empfaengerabstand in Meter (default: RMIN=1)
% NR:    Maximale Anzahl von Entfernungen (default: NR=41)
% TM:    Magnetische Dipolmoment in Ampere*Meter^2 (default:
%        TM=1)
%
% Output:
% In Abhaengigkeit von IOUT:
% FIELDS(1,NR): HZ
% FIELDS(2,NR): HR
% FIELDS(3,NR): EPHI
% HZ:    Vertikalkomponente des Magnetfeldes in Ampere/Meter,
%        Typ COMPLEX
% HR:    Radialkomponente des Magnetfeldes in Ampere/Meter,
%        Typ COMPLEX
% EPHI:  Tangentialkomponente des elektrischen Feldes in
%        Volt/Meter, Typ COMPLEX
% RR:    Aus RMIN und NRMAX gebildete Entfernungen RR(NR)
%
% Beispiel:
%  [FIELDS,RR]=VMD_F([1,0,0],1000,300,[],0,0,1,21,1);
%  berechnet HZ-Komponente fuer F=1000 Hertz, 300 Ohm*m Halbraum,
%  Sender in z=0, Empfaenger in z=0, RMIN=1 Meter, 21 Abstaende,
%  Magnetisches Dipolmoment 1 Ampere*Meter^2
%
% Aufgerufene Funktionen: HANKELFC, BTP, DOWNWARD
%
% R.-U. Boerner
%

% Defaultparameter
if(nargin<9),
    tm=1;
end
if(nargin<8),
    nr=41;
end
if(nargin<7),
    rmin=1;
end
if(nargin<6),
    zs=0;
end
if(nargin<5),
    ze=0;
end

% Ausgabeparameter:
bits=2.^(0:2)';
ifeld=sum(bits.*iout(:));
ihz=bitand(ifeld,1);
ihr=bitand(ifeld,2);
ief=bitand(ifeld,3);

if(iout(1)==1), ihz=1; end
if(iout(2)==1), ihr=1; end
if(iout(3)==1), ief=1; end

% Spez. Widerstaende, Anzahl Schichten
rho=rho(:); nl=length(rho);
% Konstanten
mu0=4e-7*pi; sigma1=1/rho(1);
% Hilfsgroessen
omega=2*pi*f;
cc=1i*omega*mu0;
c1=cc*sigma1;
% Ermittlung Sender- und Emfaengerhoehe
if(zs>0)
    warning('Sender in z>0!');
    return
end
he=ze;
hs=zs;
zm=hs-he;
zp=hs+he;

% log. aequidist. Werte fuer RR mit 10 Werten pro Dekade
q=10^0.1;
rr=rmin*q.^((1:nr)-1);

rm=sqrt(rr.^2+zm.^2);
rp=sqrt(rr.^2+zp.^2);

% Felder und Hilfsarray belegen
fields=complex(zeros(3,nr)); % Hz, Hr, Ephi
aux=complex(zeros(3,nr));

% Filterkoeffizienten
[fc0,nc,nc0]=hankelfc(3); % 1,2,3,4: sin,cos,j0,j1
[fc1,nc,nc0]=hankelfc(4);
ncnr=nc+nr-1;

% Vakuumwerte, nicht normiert
if(ze<=0),
    fields(1,:)=(3*zp.^2-rp.^2)./rp.^5; % Hz
    fields(2,:)=3*zp.*rr./rp.^5; % Hr
    fields(3,:)=rr./rp.^3; % Ephi
end

% Bereitstellung der Wellenzahlen
nu=1:ncnr;
n=nc0-nc+nu;
q=0.1*log(10);
u=exp(-(n-1)*q)/rmin;
% Admittanzen an Halbraumgrenze fuer alle Wellenzahlen u
for ii=1:ncnr,
    if(ze<=0),
        bt(ii)=btp(u(ii),f,1,rho,d);
    else
        %             bt(ii)=btp(u(ii),f,1,rho,d);
        [aa(ii),aap(ii),bt(ii)]=downward(u(ii),f,rho,d,ze);
    end
end
% Kernfunktionen
if(ze<=0),
    e=exp(u.*ze).*exp(u.*zs);
    delta=(bt-u)./(bt+u).*e;
else
    e=exp(u.*zs);
    delta=2*u.^2./(bt+u).*e;
end

% Faltung
for n=1:nr,
    for nn=1:nc,
        nu=nn+n-1;
        mn=nc-nn+1;
        nnn=nc0-nc+nu;
        u=exp(-(nnn-1)*q)/rmin;
        del0=delta(nu);
        if(ze<=0)
            del1=del0*u;
            del2=del1*u;
            aux(1,n)=aux(1,n)+del2*fc0(mn); % Hz
            aux(2,n)=aux(2,n)+del2*fc1(mn); % Hr
            aux(3,n)=aux(3,n)+del1*fc1(mn); % Ephi
        else
            del1=del0*bt(nu)*aap(nu);
            del2=del0*u*aa(nu);
            del3=del0*aa(nu);
            aux(1,n)=aux(1,n)+del2*fc0(mn); % Hz
            aux(2,n)=aux(2,n)+del1*fc1(mn); % Hr
            aux(3,n)=aux(3,n)+del3*fc1(mn); % Ephi
        end
    end
    aux(:,n)=aux(:,n)/rr(n);
end

% Absolute Feldwerte
if(ze<=0), % Lufthalbraum
    fields(1,:)=fields(1,:)-aux(1,:); % Hz
    fields(2,:)=fields(2,:)+aux(2,:); % Hr
    fields(3,:)=fields(3,:)-aux(3,:); % Ephi
else % leitfaehiger Halbraum
    fields(1,:)=aux(1,:); % Hz
    fields(2,:)=aux(2,:); % Hr
    fields(3,:)=aux(3,:); % Ephi
end

bh=1/(4*pi)+0*100; % A/m oder nT

% Normierung
fields(1,:)=tm*bh*fields(1,:); % Hz war +
fields(2,:)=tm*bh*fields(2,:); % Hr war +
fields(3,:)=-tm*cc/(4*pi)*fields(3,:); % Ephi war -

function b=btp(u,f,typ,rho,d)
% FUNCTION B=BTP(U,F,TYP,RHO,D)
% Admittance of a layered halfspace
% for TE (TYP=1) and TM (TYP=2) mode.
% U: wave number (1/m)
% F: frequency (1/s)
% RHO: layer resitivities
% D: layer thicknesses
if(nargin~=5),
    error('Wrong number of arguments.');
end
nl=length(rho);
% c=complex(0,8e-7)*pi^2*f;
mu0=4e-7*pi;
c=1i*mu0*2*pi*f;
b=sqrt(u.^2+c/rho(nl));
if(nl>1),
    beta=1;
    for nn=nl-1:-1:1
        alpha=sqrt(u.^2+c/rho(nn));
        if(typ==2),
            beta=rho(nn)/rho(nn+1);
        end
        cth=exp(-2*d(nn)*alpha);
        cth=(1-cth)./(1+cth);
        b=(b+alpha.*beta.*cth)./(beta+cth.*b./alpha);
    end
end

function [fc,nc,nc0]=hankelfc(order)
% [FC,NC,NC0]=HANKELFC(ORDER)
% Filter coefficients for Hankel transform
% FC(NC0) refers to zero argument
% NC number of coefficients
% ORDER=1: NY=+0.5  (SIN)
% ORDER=2: NY=-0.5  (COS)
% ORDER=3: NY=0.0   (J0)
% ORDER=4: NY=1.0   (J1)
% 10 data points per decade

switch(order)
    case(1)
        fc=[2.59526236E-07 3.66544843E-07 5.17830795E-07 7.31340622E-07  ...
            1.03322805E-06 1.45918500E-06 2.06161065E-06 2.91137793E-06  ...
            4.11357863E-06 5.80876420E-06 8.20798075E-06 1.15895083E-05  ...
            1.63778560E-05 2.31228459E-05 3.26800649E-05 4.61329334E-05  ...
            6.52101085E-05 9.20390575E-05 1.30122935E-04 1.83620431E-04  ...
            2.59656626E-04 3.66311982E-04 5.18141184E-04 7.30717340E-04  ...
            1.03392184E-03 1.45742714E-03 2.06292302E-03 2.90599911E-03  ...
            4.11471902E-03 5.79042763E-03 8.20004722E-03 1.15192930E-02  ...
            1.63039133E-02 2.28257757E-02 3.22249222E-02 4.47864328E-02  ...
            6.27329625E-02 8.57059100E-02 1.17418314E-01 1.53632655E-01    ...
            1.97717964E-01 2.28849849E-01 2.40311038E-01 1.65409220E-01  ...
            2.84701476E-03 -2.88016057E-01 -3.69097406E-01 -2.50107514E-02  ...
            5.71811256E-01 -3.92261572E-01 7.63280044E-02 5.16233994E-02 ...
            -6.48012082E-02 4.89047141E-02 -3.26936331E-02 2.10539842E-02  ...
            -1.33862549E-02 8.47124695E-03 -5.35123972E-03 3.37796651E-03  ...
            -2.13174466E-03 1.34513833E-03 -8.48749612E-04 5.35531006E-04  ...
            -3.37898780E-04 2.13200109E-04 -1.34520273E-04 8.48765787E-05  ...
            -5.35535069E-05 3.37899801E-05 -2.13200365E-05 1.34520337E-05  ...
            -8.48765949E-06 5.35535110E-06 -3.37899811E-06 2.13200368E-06  ...
            -1.34520338E-06 8.48765951E-07 -5.35535110E-07 3.37899811E-07];
        fc=fc(:);
        nc=80;
        nc0=40;
    case(2)
        fc=[1.63740363E-07  1.83719709E-07  2.06136904E-07  2.31289411E-07  ...
            2.59510987E-07  2.91176117E-07  3.26704977E-07  3.66569013E-07  ...
            4.11297197E-07  4.61483045E-07  5.17792493E-07  5.80972733E-07  ...
            6.51862128E-07  7.31401337E-07  8.20645798E-07  9.20779729E-07  ...
            1.03313185E-06  1.15919300E-06  1.30063594E-06  1.45933752E-06  ...
            1.63740363E-06  1.83719709E-06  2.06136904E-06  2.31289411E-06  ...
            2.59510987E-06  2.91176117E-06  3.26704977E-06  3.66569013E-06  ...
            4.11297197E-06  4.61483045E-06  5.17792493E-06  5.80972733E-06  ...
            6.51862128E-06  7.31401337E-06  8.20645798E-06  9.20779729E-06  ...
            1.03313185E-05  1.15919300E-05  1.30063594E-05  1.45933752E-05  ...
            1.63740363E-05  1.83719709E-05  2.06136904E-05  2.31289411E-05  ...
            2.59510987E-05  2.91176117E-05  3.26704977E-05  3.66569013E-05  ...
            4.11297197E-05  4.61483045E-05  5.17792493E-05  5.80972733E-05  ...
            6.51862128E-05  7.31401337E-05  8.20645798E-05  9.20779729E-05  ...
            1.03313185E-04  1.15919300E-04  1.30063594E-04  1.45933752E-04  ...
            1.63740363E-04  1.83719709E-04  2.06136904E-04  2.31289411E-04  ...
            2.59510987E-04  2.91176117E-04  3.26704976E-04  3.66569013E-04  ...
            4.11297197E-04  4.61483045E-04  5.17792493E-04  5.80972733E-04  ...
            6.51862127E-04  7.31401337E-04  8.20645797E-04  9.20779730E-04  ...
            1.03313185E-03  1.15919300E-03  1.30063593E-03  1.45933753E-03  ...
            1.63740362E-03  1.83719710E-03  2.06136901E-03  2.31289411E-03  ...
            2.59510977E-03  2.91176115E-03  3.26704948E-03  3.66569003E-03  ...
            4.11297114E-03  4.61483003E-03  5.17792252E-03  5.80972566E-03  ...
            6.51861416E-03  7.31400728E-03  8.20643673E-03  9.20777603E-03  ...
            1.03312545E-02  1.15918577E-02  1.30061650E-02  1.45931339E-02  ...
            1.63734419E-02  1.83711757E-02  2.06118614E-02  2.31263461E-02  ...
            2.59454421E-02  2.91092045E-02  3.26529302E-02  3.66298115E-02  ...
            4.10749753E-02  4.60613861E-02  5.16081994E-02  5.78193646E-02  ...
            6.46507780E-02  7.22544422E-02  8.03873578E-02  8.92661837E-02  ...
            9.80670729E-02  1.07049506E-01  1.13757572E-01  1.18327217E-01  ...
            1.13965041E-01  1.00497783E-01  6.12958082E-02 -1.61234222E-04  ...
            -1.11788551E-01 -2.27536948E-01 -3.39004453E-01 -2.25128800E-01  ...
            8.98279919E-02  5.12510388E-01 -1.31991937E-01 -3.35136479E-01  ...
            3.64868100E-01 -2.34039961E-01  1.32085237E-01 -7.56739672E-02  ...
            4.52296662E-02 -2.78297002E-02  1.73727753E-02 -1.09136894E-02  ...
            6.87397283E-03 -4.33413470E-03  2.73388730E-03 -1.72477355E-03  ...
            1.08821012E-03 -6.86602007E-04  4.33213523E-04 -2.73338487E-04  ...
            1.72464733E-04 -1.08817842E-04  6.86594042E-05 -4.33211523E-05  ...
            2.73337984E-05 -1.72464607E-05  1.08817810E-05 -6.86593962E-06  ...
            4.33211503E-06 -2.73337979E-06  1.72464606E-06 -1.08817810E-06  ...
            6.86593961E-07 -4.33211503E-07  2.73337979E-07 -1.72464606E-07];
        fc=fc(:);
        nc=164;
        nc0=122;
    case(3)
        fc=[2.89878288E-07  3.64935144E-07  4.59426126E-07  5.78383226E-07 ...
            7.28141338E-07  9.16675639E-07  1.15402625E-06  1.45283298E-06  ...
            1.82900834E-06  2.30258511E-06  2.89878286E-06  3.64935148E-06  ...
            4.59426119E-06  5.78383236E-06  7.28141322E-06  9.16675664E-06  ...
            1.15402621E-05  1.45283305E-05  1.82900824E-05  2.30258527E-05  ...
            2.89878259E-05  3.64935186E-05  4.59426051E-05  5.78383329E-05  ...
            7.28141144E-05  9.16675882E-05  1.15402573E-04  1.45283354E-04  ...
            1.82900694E-04  2.30258630E-04  2.89877891E-04  3.64935362E-04  ...
            4.59424960E-04  5.78383437E-04  7.28137738E-04  9.16674828E-04  ...
            1.15401453E-03  1.45282561E-03  1.82896826E-03  2.30254535E-03  ...
            2.89863979E-03  3.64916703E-03  4.59373308E-03  5.78303238E-03  ...
            7.27941497E-03  9.16340705E-03  1.15325691E-02  1.45145832E-02  ...
            1.82601199E-02  2.29701042E-02  2.88702619E-02  3.62691810E-02  ...
            4.54794031E-02  5.69408192E-02  7.09873072E-02  8.80995426E-02  ...
            1.08223889E-01  1.31250483E-01  1.55055715E-01  1.76371506E-01  ...
            1.85627738E-01  1.69778044E-01  1.03405245E-01 -3.02583233E-02  ...
            -2.27574393E-01 -3.62173217E-01 -2.05500446E-01  3.37394873E-01  ...
            3.17689897E-01 -5.13762160E-01  3.09130264E-01 -1.26757592E-01  ...
            4.61967890E-02 -1.80968674E-02  8.35426050E-03 -4.47368304E-03  ...
            2.61974783E-03 -1.60171357E-03  9.97717882E-04 -6.26275815E-04  ...
            3.94338818E-04 -2.48606354E-04  1.56808604E-04 -9.89266288E-05  ...
            6.24152398E-05 -3.93805393E-05  2.48472358E-05 -1.56774945E-05  ...
            9.89181741E-06 -6.24131160E-06  3.93800058E-06 -2.48471018E-06  ...
            1.56774609E-06 -9.89180896E-07  6.24130948E-07 -3.93800005E-07  ...
            2.48471005E-07 -1.56774605E-07  9.89180888E-08 -6.24130946E-08];
        fc=fc(:);
        nc=100;
        nc0=60;
    case(4)
        fc=[1.84909557E-13  2.85321327E-13  4.64471808E-13  7.16694771E-13 ...
            1.16670043E-12  1.80025587E-12  2.93061898E-12  4.52203829E-12 ...
            7.36138206E-12  1.13588466E-11  1.84909557E-11  2.85321327E-11 ...
            4.64471808E-11  7.16694771E-11  1.16670043E-10  1.80025587E-10 ...
            2.93061898E-10  4.52203829E-10  7.36138206E-10  1.13588466E-09 ...
            1.84909557E-09  2.85321326E-09  4.64471806E-09  7.16694765E-09 ...
            1.16670042E-08  1.80025583E-08  2.93061889E-08  4.52203807E-08 ...
            7.36138149E-08  1.13588452E-07  1.84909521E-07  2.85321237E-07 ...
            4.64471580E-07  7.16694198E-07  1.16669899E-06  1.80025226E-06 ...
            2.93060990E-06  4.52201549E-06  7.36132477E-06  1.13587027E-05 ...
            1.84905942E-05  2.85312247E-05  4.64449000E-05  7.16637480E-05 ...
            1.16655653E-04  1.79989440E-04  2.92971106E-04  4.51975783E-04 ...
            7.35565435E-04  1.13444615E-03  1.84548306E-03  2.84414257E-03 ...
            4.62194743E-03  7.10980590E-03  1.15236911E-02  1.76434485E-02 ...
            2.84076233E-02  4.29770596E-02  6.80332569E-02  9.97845929E-02 ...
            1.51070544E-01  2.03540581E-01  2.71235377E-01  2.76073871E-01 ...
            2.16691977E-01 -7.83723737E-02 -3.40675627E-01 -3.60693673E-01 ...
            5.13024526E-01 -5.94724729E-02 -1.95117123E-01  1.99235600E-01 ...
            -1.38521553E-01  8.79320859E-02 -5.50697146E-02  3.45637848E-02 ...
            -2.17527180E-02  1.37100291E-02 -8.64656417E-03  5.45462758E-03 ...
            -3.44138864E-03  2.17130686E-03 -1.36998628E-03  8.64398952E-04 ...
            -5.45397874E-04  3.44122545E-04 -2.17126585E-04  1.36997597E-04 ...
            -8.64396364E-05  5.45397224E-05 -3.44122382E-05  2.17126544E-05 ...
            -1.36997587E-05  8.64396338E-06 -5.45397218E-06  3.44122380E-06 ...
            -2.17126543E-06  1.36997587E-06 -8.64396337E-07  5.45397218E-07];
        fc=fc(:);
        nc=100;
        nc0=60;
    otherwise
        error('Order must be 1<=order<=4.');
end

