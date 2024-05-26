function [InvSqrtMat] = InvSqrtMatrix(MM); 

[eigvec,eigval] = eigs(MM);

SqrtEigval = diag((diag(abs(eigval))).^.5);
InvSqrtMat = eigvec*SqrtEigval;

