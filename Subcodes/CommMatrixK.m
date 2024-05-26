function [mK]=CommMatrixK(k)

	%mK = unit(k*k)[vec(reshape(range(1,k*k)-1,k,k))][];
    
    range=[1:k*k];
    sel=vec(reshape(range,k,k)');
    II=eye(k*k);
    mK = II(sel,:);

    
end