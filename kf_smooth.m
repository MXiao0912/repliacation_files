function [alfa, var] = kf_smooth(v, invF, L, a, P, Z)

TT = size(v,2);
m = size(P,1);

alfa = NaN(m,TT);
var = NaN(m,m,TT);
r = NaN(m, TT);
N = NaN(m,m,TT);

r(:,TT) = 0;
N(:,:,TT) = 0;

for i=TT:-1:2
    r(:,i-1) = Z{i}'* invF{i}*v{i}+L(:,:,i)'*r(:,i);
    N(:,:,i-1) = Z{i}'*invF{i}*Z{i} + L(:,:,i)'*N(:,:,i)*L(:,:,i);
    alfa(:,i) = a(:,i)+P(:,:,i)*r(:,i-1);
    var(:,:,i) = P(:,:,i)-P(:,:,i)*N(:,:,i-1)*P(:,:,i);
end

r0 = Z{i}'* invF{1}*v{1}+L(:,:,1)'*r(:,1);
N0 = Z{i}'*invF{1}*Z{i} + L(:,:,1)'*N(:,:,1)*L(:,:,1);
alfa(:,1) = a(:,1)+P(:,:,1)*r0;
var(:,:,1) = P(:,:,1)-P(:,:,1)*N0*P(:,:,1);

end