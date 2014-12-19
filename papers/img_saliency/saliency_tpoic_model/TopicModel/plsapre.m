function [pz pdz pwz pzdw]=plsa(x,k)
% function [pz pdz pwz pzdw]=plsa(x,k)
% Return pLSA probability matrix p of m*n matrix x
% x is the document-word co-occurence matrix
% k is the number of the topics--z
% document--collums,word--rows

err=0.000001
[m n]=size(x);
pz=rand(1,k);
pz2=pz;
pz2(1)=pz2(1)+2*err;
pdz=rand(k,n);
pwz=rand(k,m);
pzdw=rand(m,n,k);        %initialize
h=0;

deno=zeros(1,k);        %denominator of p(d/z) and p(w/z)
denopzdw=zeros(m,n);       %denominator of p(z/d,w)
numepdz=zeros(k,n);        %numerator of p(d/z)
numepwz=zeros(k,m);        % numerator of p(w/z)
R=sum(sum(x))

for ki=1:k
for i=1:m
    for j=1:n
        deno(ki)=deno(ki)+x(i,j)*pzdw(i,j,ki);
    end
end
end
% p(d/z)
for ki=1:k
    for j=1:n
        for i=1:m
            numepdz(ki,j)=numepdz(ki,j)+x(i,j)*pzdw(i,j,ki);
        end
        pdz(ki,j)=numepdz(ki,j)/deno(ki);
    end
end
% disp(pdz);

% p(w/z)
for ki=1:k
    for i=1:m
        for j=1:n
            numepwz(ki,i)=numepwz(ki,i)+x(i,j)*pzdw(i,j,ki);
        end
        pwz(ki,i)=numepwz(ki,i)/deno(ki);
    end
end
% disp(pwz);

% p(z)
for ki=1:k
    pz(ki)=deno(ki)./R;
end
%denominator of p(z/d,w)
for i=1:m
    for j=1:n
        for ki=1:k
            denopzdw(i,j)=denopzdw(i,j)+pz(ki)*pdz(ki,j)*pwz(ki,i);
        end
    end
end
% p(z/d,w)
for i=1:m
     for j=1:n
         for ki=1:k
             pzdw(i,j,ki)=pz(ki)*pdz(ki,j)*pwz(ki,i)/denopzdw(i,j);
         end
     end
end
% fprintf('p(z/d,w)=\n');
% disp(pzdw)

%  iteration 
iteration=0;
fprintf('iteration:\n');
while abs(pz2(1)-pz(1))>err | abs(pz2(2)-pz(2))>err
    iteration=iteration+1;
    deno=zeros(1,k);        %denominator of p(d/z) and p(w/z)
    denopzdw=zeros(m,n);       %denominator of p(z/d,w)
    numepdz=zeros(k,n);        %numerator of p(d/z)
    numepwz=zeros(k,m);        % numerator of p(w/z)
    fprintf('iteration %d:\n',iteration);
    for ki=1:k
for i=1:m
    for j=1:n
        deno(ki)=deno(ki)+x(i,j)*pzdw(i,j,ki);
    end
end
end
% p(d/z)
for ki=1:k
    for j=1:n
        for i=1:m
            numepdz(ki,j)=numepdz(ki,j)+x(i,j)*pzdw(i,j,ki);
        end
        pdz(ki,j)=numepdz(ki,j)/deno(ki);
    end
end
fprintf('p(d/z)=\n');
disp(pdz)
% p(w/z)
for ki=1:k
    for i=1:m
        for j=1:n
            numepwz(ki,i)=numepwz(ki,i)+x(i,j)*pzdw(i,j,ki);
        end
        pwz(ki,i)=numepwz(ki,i)/deno(ki);
    end
end
fprintf('p(w/z)=\n');
disp(pwz)

% p(z)
pz=pz2;
for ki=1:k
    pz2(ki)=deno(ki)./R;
end
fprintf('p(z)=\n');
disp(pz2)
%denominator of p(z/d,w)
for i=1:m
    for j=1:n
        for ki=1:k
            denopzdw(i,j)=denopzdw(i,j)+pz2(ki)*pdz(ki,j)*pwz(ki,i);
        end
    end
end
% p(z/d,w)
for i=1:m
     for j=1:n
         for ki=1:k
             pzdw(i,j,ki)=pz2(ki)*pdz(ki,j)*pwz(ki,i)/(denopzdw(i,j)+eps);
         end
     end
end
end     %end while
return;




