function [x]=fullProj(y0,k,n,x)
% taken from:
% https://github.com/canyilu/Projection-onto-the-capped-simplex
%
% Inputs
%   y0: nx1 vector.
%   k: the sum constraint.
%   n: dimension of vector
%   x: all zero vector for faster indexing
%
% Outputs
%   x: Nx1 vector.
%
% Reference:
%   Weiran Wang: "Projection onto the capped simplex".
%     March 3, 2015, arXiv:1503.01002.
%
% Weiran Wang. 03/03/2015.

% x = zeros(q,1);
if k==n
    x=ones(n,1);
    return;
end
[y,idx]=sort(y0,'ascend');

% Test the possiblity of a==b are integers.
if k==round(k)
    b=n-k;
    if y(b+1)-y(b)>=1
        x(idx(b+1:end))=1;
        return;
    end
end

% Assume a=0.
s=cumsum(y);
y=[y;inf];
for b=1:n
    % Hypothesized gamma.
    gamma = (k+b-n-s(b)) / b;
    if ((y(1)+gamma)>0) && ((y(b)+gamma)<1) && ((y(b+1)+gamma)>=1)
        x(idx)=[y(1:b)+gamma; ones(n-b,1)];
        return;
    end
end

% Now a>=1;
for a=1:n
    for b=a+1:n
        % Hypothesized gamma.
        gamma = (k+b-n+s(a)-s(b))/(b-a);
        if ((y(a)+gamma)<=0) && ((y(a+1)+gamma)>0) && ((y(b)+gamma)<1) && ((y(b+1)+gamma)>=1)
            x(idx)=[zeros(a,1); y(a+1:b)+gamma; ones(n-b,1)];
            return;
        end
    end
end

% edge case: 0 < abs( k - round(k) ) < epsilon
xtmp = [zeros(n-ceil(k),1); k-floor(k); ones(floor(k),1)];
x(idx) = xtmp;
end
