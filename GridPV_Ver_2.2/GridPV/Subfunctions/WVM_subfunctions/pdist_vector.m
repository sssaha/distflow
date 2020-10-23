function D = pdist_vector(X)

sz=size(X);
if sz(2)==1
    X=X';
end



% D = squareform(sqrt(bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y)));

%elminate need for squareform function (and, hence, Statistics Toolbox)
D1=sqrt(bsxfun(@plus,dot(X,X,1)',dot(X,X,1))-2*(X'*X));
T1=tril(ones(size(D1)),-1)~=0;
D=D1(T1);
D=D(:)';

%Adapted from http://statinfer.wordpress.com/2011/11/14/efficient-matlab-i-pairwise-distances/