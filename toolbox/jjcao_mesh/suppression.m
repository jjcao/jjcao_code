function newFun = suppression(F,oldFunction)

%The suppression operator used in the computation of the mesh saliency based on the ACM SIGGRAPH 2005
%paper "Mesh Saliency"

%Hui Wang, June, 5, 2011, wanghui19841109@gmail.com

MAX = max(oldFunction);
MIN = min(oldFunction);

newFun = (oldFunction - MIN) / (MAX - MIN);
[minimaIndex,maximaIndex] = extremumIndex(F,newFun);
value = newFun(maximaIndex);
average = (sum(value) - 1.0) / (length(value) - 1);
factor = (1.0 - average) .^2;
newFun = factor * newFun;

