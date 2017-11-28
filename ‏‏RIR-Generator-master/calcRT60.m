function [RT60,V] = calcRT60(roomDim,beta)
V = prod(abs(roomDim));
S(1:2) =  roomDim(1)*roomDim(2);
S(3:4) =  roomDim(1)*roomDim(3);
S(5:6) =  roomDim(2)*roomDim(3);
if length(beta)==1
    beta = beta*ones(6,1);
else
    beta = beta(:);
end
c = 340;
numer = 24*log(10)*V;
denom = c*S*(1-beta.^2);
RT60 = numer/denom;

end