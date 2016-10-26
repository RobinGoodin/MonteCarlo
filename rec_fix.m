function y = rec_fix(n, k, s)
if(n == 1)
    y = fix(k / s);
else
    y = rec_fix(n - 1, fix(k / s) , s);
end