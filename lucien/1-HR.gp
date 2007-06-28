default(parisize, 100M)
default(realprecision, 76)

global(A,B,X0,X1,Y0,Y1,x0,x1,y0,y1)

Y0=4*X0*X1*(X1-A*X0)*(X1+B*X0);
Y1=(X1^2+A*B*X0^2)^2;

doe(xx)=
{
[subst(subst(y0,x0,xx[1]),x1,xx[2]),
subst(subst(y1,x0,xx[1]),x1,xx[2])]
}

Around_the_top(aa,bb,k,n,cutoff)=
{
local(einde,xx,uit);
einde=0;
xx=vector(2);
if(n,,einde=1);
if(k,,einde=1);
if(einde,print("Usage: Around_the_horn(a,b,k,n)\n"
"where a,b are from abc,\n"
"k is the number of iderations, and \n"
"n is the number of steps.");return());
y0=subst(subst(subst(subst(Y0,A,aa),B,bb),X0,x0),X1,x1);
y1=subst(subst(subst(subst(Y1,A,aa),B,bb),X0,x0),X1,x1);
for(i=0,n,\
	angle=i*Pi/n;\
	xx=[cos(angle),sin(angle)];\
	for(j=1,k,xx=doe(xx));\
	uit=log(sqrt(xx[1]^2+xx[2]^2))/4^k;\
	if(uit<cutoff,print(i));\
	print(uit);\
);
return()
}

