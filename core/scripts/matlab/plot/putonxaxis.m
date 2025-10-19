function putonxaxis(pos,mark,mindist,takey)

% PUTONXAXIS - Puts string on x axis or replaces it
% putonxaxis(position,string) puts string on position pos
% putonxaxis(position) includes tick at position
% putonxaxis(position,string,mindist) deletes ticks close (mindist) to position

if nargin<1, error('Specify position!'); end
if nargin<2, mark=num2str(pos); end
if nargin<3, mindist=0; end
if nargin<4, takey=0; end

if takey,
    xt=get(gca,'YTick');
    xtl=get(gca,'YTickLabel');
else
    xt=get(gca,'XTick');
    xtl=get(gca,'XTickLabel');
end
l=1;
newxt=[];newxtl={};
for i=1:length(xt),
    if pos==xt(i), %replace
        newxtl{l}=mark;
        newxt(l)=pos;
    else
        if abs(pos-xt(i))>mindist, % take it
            newxt(l)=xt(i);
            if ischar(xtl), newxtl{l}=xtl(i,:); else newxtl{l}=xtl{i}; end        
        else 
            l=l-1;
        end
        if (i<length(xt))&&(pos>xt(i))&(pos<xt(i+1)), % insert
            l=l+1;
            newxt(l)=pos;
            newxtl{l}=mark;
        end
    end
    l=l+1;    
end
if (pos>xt(end))&(abs(pos-xt(i))>mindist), newxt(end+1)=pos;newxtl{end+1}=mark; end
if takey,
    set(gca,'XTick',newxt,'XTickLabel',newxtl);
else
    set(gca,'XTick',newxt,'XTickLabel',newxtl);
end
