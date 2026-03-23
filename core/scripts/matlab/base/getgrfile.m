function [outfile,iseps,ispdf]=getgrfile(infile)

% GETGRFILE - get graphics file
% calls UIGETFILE to determine file name and type (eps or png)
% outfile = getgrfile(basefilename);
% [outfile,iseps] = getgrfile(basefilename);

[path,name,ext]=fileparts(infile); 
outfile=strrep(infile,ext,'.png');
[fname,pname]=uiputfile({'*.png ';'*.eps ';'*.pdf '},'Save Figure as',outfile);
ispdf=any(strfind(fname,'.pdf'));
iseps=any(strfind(fname,'.eps'))|ispdf;
if fname==0,
    outfile='';
else
    if ispdf, 
        ext='.eps'; 
        fname=strrep(fname,'.pdf','.eps');
    end
    outfile=fullfile(pname,fname);
    [ff,pp,ee]=fileparts(outfile);
    if strcmp(ee,''), 
        ext='.png';
        if iseps, ext='.eps'; end
        outfile=[outfile ext]; 
    end
end
