function ALLTEM = readusffile(filename,nr)

% READUSFFILE - Read universal sounding file (usf)
% TEM = readusffile(filename)

if nargin<2, nr=0; end
fid=fopen(filename,'r');
zeile=fgetl(fid);
if strcmp(zeile(1:2),'//'),
    dpos=strfind(zeile,':');
    tok=upper(zeile(3:dpos-1));
    while ~strcmp(tok,'SOUNDINGS'),
        zeile=fgetl(fid);
        dpos=strfind(zeile,':');
        tok=upper(zeile(3:dpos-1));
    end
    ns=str2num(zeile(dpos+1:end));
end
for i=1:ns,
    TEM=[];TEM.ndata=0;tok='';
    while isstr(zeile)&&(~isequal(zeile,'/END')),
        if (length(zeile)>0)&&(zeile(1)=='/'),
            dpos=strfind(zeile,':');
            if isnumeric(dpos),
                tok=upper(zeile(2:dpos-1));
                if isequal(tok,'POINTS'), TEM.ndata=str2num(zeile(dpos+1:end)); end
                if isequal(tok,'CURRENT'), TEM.current=str2num(zeile(dpos+1:end)); end
                if isequal(tok,'COIL_SIZE'), TEM.coilsize=str2num(zeile(dpos+1:end)); end
                if isequal(tok,'LOOP_SIZE'), TEM.loopsize=str2num(strrep(zeile(dpos+1:end),',',' ')); end
                if isequal(tok,'RAMP_TIME'), TEM.ramptime=str2num(strrep(zeile(dpos+1:end),',',' ')); end
                if isequal(tok,'TIME_DELAY'), TEM.timedelay=str2num(strrep(zeile(dpos+1:end),',',' ')); end
            end
        end
        zeile=fgetl(fid);
    end
    if TEM.ndata>0,
        fstr='%d,%f,%f';
        it=2;iv=3;ie=4;
        zeile=strrep(fgetl(fid),',','');
        %% tokenzeile auswerten
        zeile=strrep('INDEX,	TIME,	VOLTAGE	ST_DEV',',','');
        [tok,zeile]=strtok(zeile);
        fstr='';
        rows=0;
        while length(tok)>0,
            fstr=[fstr '%f,'];
            rows=rows+1;
            [tok,zeile]=strtok(zeile);
        end
        fstr(end)='';
        %%
        A=fscanf(fid,fstr,[rows,TEM.ndata])';
        if it, TEM.t=A(:,it); end
        if iv, TEM.v=A(:,iv); end
        if ie<=rows, TEM.e=A(:,ie); end
    end
    while ~isequal(zeile,'/END'), 
        zeile=fgetl(fid); 
    end
    zeile=fgetl(fid);
    if nr>0,
        if i==nr, ALLTEM=TEM;break; end
    else
        ALLTEM{i}=TEM;
    end
end
fclose(fid);
