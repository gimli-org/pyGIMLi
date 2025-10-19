function [stringline,comment]=destrip(stringline,commentchar)

% DESTRIP - STRIPS comment from string
% remainder = destrip(string)
% [remainder,comment] = destrip(string)
% ... = destrip(string,commentcharacter)

if nargin<2, commentchar='#'; end
aa=strfind(stringline,commentchar);
if ~isempty(aa), stringline=stringline(1:aa(1)-1); end
if nargout>1, comment=stringline(aa(1)+1:end); end
