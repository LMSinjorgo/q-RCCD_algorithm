function [x]=fullProj(y0,k,n,x)
% function [x,e]=projection(y0,k) solves the projection onto
%   the capped simplex:
%     min_x 0.5*|x-y0|^2
%       s.t. 0 <= x_i <= 1,
%            sum_i x_i = k.
%
% Inputs
%   y0: nx1 vector.
%   k: the sum constraint.
%
% Outputs
%   x: Nx1 vector.
%   e: 0.5*|x-y0|^2.
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
xtmp = [zeros(n-floor(k)-1,1); k-floor(k); ones(floor(k),1)];
x(idx) = xtmp;
end