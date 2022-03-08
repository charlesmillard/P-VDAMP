function out = nmse(x,x0)
    x_mnsq = mean(abs(x0(:)).^2);
    out = 10*log10(mean(abs(x(:) - x0(:)).^2)./x_mnsq);
end

