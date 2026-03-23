function chiq=chi2(soll,ist,err,lolo)

% CHI2 - computation of chi^2-error
% chiq = chi2(data,response,error[,islog])
% data, response .. vectors to compare
% error .. error/standard deviation
% islog .. treat data and response as logarithmic (default=1)

if nargin<3, err=0.03; end
if nargin<4, lolo=1; end
if lolo,
%     chiq=mean(((log(ist)-log(soll))./log(1+err)./log(soll)).^2);
    chiq=mean(((log(ist(:))-log(soll(:)))./log(1+err(:))).^2);
else
%     chiq=mean(((ist(:)./soll(:)-1).^2)./(err(:).^2));
    chiq=mean(((ist(:)-soll(:))./err(:)).^2);
end