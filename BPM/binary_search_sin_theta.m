function x = binary_search_sin_theta(n1, n2, d, lambda)

f = @(x)(n1 / n2 * x - sin(pi * d / lambda * x));

a = 0.0;
b = lambda / d;

while(abs(a-b) > 1e-14)
    x = 0.5 * (a + b);
    if (f(x)>0)
        b = x;
    elseif (f(x)<0)
        a = x;
    else
        break;
    end
end

end

