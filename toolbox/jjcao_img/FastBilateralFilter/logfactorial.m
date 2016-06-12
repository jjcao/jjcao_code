% aproximate log-factorial using Sterling's formula
function val = logfactorial(n)
   val = n*(log(n)-1) + 0.5*log(2*pi*n);
end