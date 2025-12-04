function y = f( a, x )
% equivalent to std. tanh(), but faster

y = 2 ./ ( 1 + exp( x .* (-2) ) ) - 1;