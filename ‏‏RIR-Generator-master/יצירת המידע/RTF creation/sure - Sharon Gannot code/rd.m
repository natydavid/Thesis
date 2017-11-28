function [x,f] = rd(name,machine)

%[x,f] = rd('bla','l')

if nargin == 1 machine = 'native'; end;
f = fopen(name,'r',machine);
if f ~= -1 
  x = fread(f,'short'); fclose(f); 
else 
  warning('No such file'); 
  x = nan;
end
