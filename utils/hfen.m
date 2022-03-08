function h = hfen(a, a0)
    % compute High-Frequency Error Norm (HFEN)
    f = fspecial('log',15,1.5);
    
    blur = conv2(a-a0, f);
    
    h = 10*log10(norm(blur(:)).^2/norm(a0(:)).^2);
end

