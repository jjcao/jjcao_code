function  T  =   maxFilter(f0 , w)
%
%  compute the maximum local dynamic range
%
%  T = max_x  max_{|y| <= r} |f(x) - f(x-y)|
%
%  f0     : grayscale image
%  w           : search window (w = 2*r+1, must be odd)
%  T           : maximum local dynamic range
%
T = -1;
sym    = (w - 1)/2;
[m, n] = size(f0);
template = f0;

% scan along row
for ii = 1 : m
    L     = zeros(n, 1);
    R     = zeros(n, 1);
    L(1)  = template(ii, 1);
    R(n)  = template(ii, n);
    
    for k = 2 : n
        if  mod(k - 1, w) == 0
            L(k)          = template(ii ,  k);
            R(n - k + 1)  = template(ii ,  n - k + 1);
        else
            L(k)          = max( L(k-1) , template(ii, k) );
            R(n - k + 1)  = max( R(n - k + 2), template(ii, n - k + 1) );
        end
    end
    
    for k = 1 : n
        p = k - sym;
        q = k + sym;
        if p < 1
            r = -1;
        else
            r = R(p);
        end
        if q > n
            l = -1;
        else
            l = L(q);
        end
        template(ii, k) = max(r,l);
    end
    
end

% scan along column
for jj = 1 : n
    L    = zeros(m, 1);
    R    = zeros(m, 1);
    L(1) = template(1, jj);
    R(m) = template(m, jj);
    
    for k = 2 : m
        if  mod(k - 1, w) == 0
            L(k)          = template(k, jj);
            R(m - k + 1)  = template(m - k + 1, jj);
        else
            L(k)          = max( L(k - 1), template(k, jj) );
            R(m - k + 1)  = max( R(m - k + 2), template(m - k + 1, jj));
        end
    end
    
    for k = 1 : m
        p = k - sym;
        q = k + sym;
        if p < 1
            r = -1;
        else
            r = R(p);
        end
        if q > m
            l = -1;
        else
            l = L(q);
        end
        temp = max(r,l) - f0(k, jj);
        if temp > T
            T = temp;
        end
    end
end
end

