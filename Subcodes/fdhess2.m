

function h = fdhess2( f, x0 )
    %local k, hessian, grdd, ax0, dax0, dh, xdh, ee, f0, i, j, eps;
    %local f:proc;

    %/* check for complex input */
    if imag(x0) ~= 0;
        if max(imag(x0)) > eps
            error('ERROR: Not implemented for complex matrices.')
        else
            x0 = real(x0);
        end
    end

%/* initializations */
    k = length(x0);
    hessian = zeros(k,k);
    grdd = zeros(k,1);
    epsil = 6.0554544523933429e-6;

%/* Computation of stepsize (dh) */
    ax0 = abs(x0);
    if x0 ~= 0;
        dax0 = x0./ax0;
    else
        dax0 = 1;
    end
    dh  = epsil*max([ax0 (1e-2)*ones(k,1)],[],2).*dax0;

    xdh = x0+dh;
    dh  = xdh-x0;    %/* This increases precision slightly */

    ee  = diag(dh);

%/* Computation of f0=f(x0) */
    f0 = feval(f,x0);

%/* Compute forward step */
    for ii=1:k;
        grdd(ii) = feval(f,x0+ee(:,ii));
    end     
        
%/* Compute "double" forward step */
  
    for ii=1:k;
        for jj=ii:k;
            hessian(ii,jj) = feval(f, x0+ee(:,ii)+ee(:,jj) );
            if ii ~= jj;
                hessian(jj,ii) = hessian(ii,jj);
            end
        end
    end
    h =( ( ( (hessian - repmat(grdd,1,k)) - repmat(grdd',k,1)) + f0 ) ./ ( dh*dh') );
end