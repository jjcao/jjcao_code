function pzw = getpzw(pwz,pz)
%m topic
%n words
[m,n] = size(pwz);
pzw = zeros(m,n);
for k = 1:m
    for j = 1:n
        pzw(k,j) = pwz(k,j)*pz(k);
        tmp = 0;
        for num = 1:m
            tmp = tmp+pwz(num,j)*pz(num);
        end
        pzw(k,j) = pzw(k,j)/tmp;
    end
end
    