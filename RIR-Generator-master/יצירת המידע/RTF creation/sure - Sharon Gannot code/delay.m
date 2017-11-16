function x_delayed = delay(x,D);
%
% Usage:  x_delayed = delay(x,D);
%
% x - the original signal
% D - amount of delay;
% if D > 0 there are D zeros in the begining and the last D samples are lost.
% if D < 0 the first D samples are lost and there are D zeros in the end.
%
% x_delayed - is the delayed version of x;

tr = 0;
[m,n] = size(x);
if m > n
  x = x';
  tr = 1;
end;

x_delayed = zeros(size(x));
	
if D > 0
  
  x_delayed(:,D+1:length(x)) = x(:,1:length(x)-D);

elseif D <= 0
  
  x_delayed(:,1:length(x)-abs(D)) = x(:,abs(D)+1:length(x));

end;

if tr, x_delayed = x_delayed'; end;