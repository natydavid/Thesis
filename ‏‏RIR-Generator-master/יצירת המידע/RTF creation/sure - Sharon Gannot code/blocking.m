function [SNRin,SNRknw,SNRnst,h_known] = blocking(s,v,h,FLE_l,FLE_r,T,Tblk);

h_known = known(s,T,FLE_l,FLE_r);

% --------------------------
% Check Blocking of the filters
% ----------------------------

[SNRin,SNRknw,SNRnst] = block(s,v,h,h_known,FLE_l,Tblk)

% ------------------------------
% Find filters from clean signal
% ------------------------------

function h_known = known(S,T,FLE_l,FLE_r);

M = size(S,1);
h_known = zeros(M,FLE_l+FLE_r);
h_known(1,FLE_l+1) = 1;

for m = 2:M
  h_known(m,:) = est2bat(S(1,T)',S(m,T)',FLE_l,FLE_r);
end;

% --------------------------
% Find blocking capability
% --------------------------

function [SNRin,SNRknw,SNRnst] = block(S,V,h,h_known,FLE_l,Tblk);

M = size(S,1);

for m = 1:M-1
  
  s_knw = S(m+1,:) - delay(filter(h_known(m+1,:),1,S(1,:)),-FLE_l); 
  v_knw = V(m+1,:) - delay(filter(h_known(m+1,:),1,V(1,:)),-FLE_l);   
  s_nst = S(m+1,:) - delay(filter(h(m+1,:),1,S(1,:)),-FLE_l); 
  v_nst = V(m+1,:) - delay(filter(h(m+1,:),1,V(1,:)),-FLE_l);   
  
  SNRknw(m) = 10*log10(var(s_knw(Tblk))/var(v_knw(Tblk)));
  SNRnst(m) = 10*log10(var(s_nst(Tblk))/var(v_nst(Tblk)));
  
end;

for m = 1:M
  
  SNRin(m)  = 10*log10(var(S(m,Tblk))/var(V(m,Tblk)));  

end;
