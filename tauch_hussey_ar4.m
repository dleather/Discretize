n	= 4;
mx	= 0;
rx	= 0.90;
s	= 0.01;

[xx,wx]=gauss_herm(n);		% nodes and weights for x
st=sqrt(2)*s*xx+mx
x=xx(:,ones(n,1));
y=x';
w=wx(:,ones(n,1))';

p=(exp(y.*y-(y-rx*x).*(y-rx*x)).*w)./sqrt(pi);
sm=sum(p')';
p=p./sm(:,ones(n,1))