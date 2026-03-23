function showmaxmindata(data,resp,freq,npl)

if nargin<4, npl=2; end
if nargin<3, freq=[110,220,440,880,1760,3520,7040,14080,28160,56320]; end


subplot(1,npl,npl-1);
semilogy(data(:,1),freq,'ro-',resp(:,1),freq,'bx-');
axis tight;grid on;title('real part');
subplot(1,npl,npl);
semilogy(data(:,2),freq,'ro-',resp(:,2),freq,'bx-');
axis tight;grid on;title('imaginary part');
