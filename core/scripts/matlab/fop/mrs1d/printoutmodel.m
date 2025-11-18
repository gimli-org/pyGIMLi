function printoutmodel(thk,wc,t2,misfit,iter)

if nargin<5, iter=0; end
fprintf('Iter: %d\n',iter);
% fprintf('wc =');fprintf('%.2f ',wc);fprintf('\n');
% fprintf('t2 =');fprintf('%.2f ',t2);fprintf('\n');
fprintf('wc =');fprintf('%.1f ',wc*100);fprintf('%%\n');
fprintf('t2 =');fprintf('%.1f ',t2*1000);fprintf('ms\n');
fprintf('thk=');fprintf('%.2f ',thk);fprintf('m\n');
fprintf('z=0 ');fprintf('%.2f ',cumsum(thk));fprintf('m\n');
rrms=norm(misfit)/sqrt(prod(size(misfit)));
fprintf('rms =%.2fnV\n',rrms);
