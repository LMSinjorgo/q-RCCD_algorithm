function [e] = quickNormEst(A, q)
% quick custom implementation of normEst function, to bypass slow matlab
% built in function.
% works only if A is nonnegative, symmetric and positive diagonal

if nnz(A) == q
   % then matrix is diagonal
   e = full(max(max(A)));
   return;
end
tol = 1.e-6;
x = full(sum(A,1))';
cnt = 0;
e = norm(x);
x = x/e;
e0 = 0;
while abs(e-e0) > tol*e
   e0 = e;
   Ax = A*x;
   if nnz(Ax) == 0
      Ax = rand(size(Ax),class(Ax));
   end
   x = A*Ax;
   normx = norm(x);
   e = normx/norm(Ax);
   x = x/normx;
   cnt = cnt+1;
   if cnt > 100
       % no clear convergence, but no matter
       break;
   end
end
end

