function epsprint(fig,outfile,topdf,bw)

% EPSPRINT - print figure into eps graphics file
% epsprint(figure_handle,filename)
% epsprint(figure_handle,filename,convert_to_pdf)
% epsprint(filename[,convert_to_pdf]) takes current figure (gcf)

if (nargin==1)&&ischar(fig), outfile=fig;fig=gcf; end
set(fig,'PaperPositionMode','manual');
set(fig,'Units',get(fig,'PaperUnits'));
po=get(fig,'PaperPosition');
newpo=po([4 3])-po([2 1]);
if isequal(get(fig,'PaperPositionMode'),'manual')&&(newpo(1)>0)&(newpo(2)>0),
    set(fig,'PaperSize',po([4 3])-po([2 1]));
    set(fig,'PaperPositionMode','auto');
end
if nargin>3, %greyscale
    print(fig,'-deps','-painters',outfile);
else
    print(fig,'-depsc2','-painters',outfile);
end
if (nargin>2)&&(topdf),
	dos(['epstopdf "' strrep(outfile,'.eps',''),'.eps"']);
end
