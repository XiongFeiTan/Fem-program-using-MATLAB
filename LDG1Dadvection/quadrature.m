function quad = quadrature(f,a,b)
%Define variables:
%   a,b,f: integrand
%   w:weight
%   A:integration point
%   quad:result

A=[-0.9602898565 -0.7966664774 -0.5255324099 -0.1834346425 ...
    0.1834346425 0.5255324099 0.7966664774 0.9602898565];

w=[0.1012285363 0.2223810345 0.3137066459 0.3626837834 ...
   0.3626837834 0.3137066459 0.2223810345 0.1012285363];

%mapping
T = ((a+b)/2) + ((b-a)/2)*A;
%calculate
quad = ((b-a)/2)*sum(w.*feval(f,T));
end

