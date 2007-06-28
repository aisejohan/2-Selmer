default(parisize, 100M)
default(realprecision, 100)

a = 125010010110000000000000000;
print("The value of a is: ",a)
b = 121230120000000912323123452;
print("The value of b is: ",b)

Y0=4*X0*X1*(X1-a*X0)*(X1+b*X0);
Y1=(X1^2+a*b*X0^2)^2;

UIT0=X0;
UIT1=X1;

compose()=
{
UIT0=subst(subst(UIT0,X0,XX0),X1,XX1);
UIT1=subst(subst(UIT1,X0,XX0),X1,XX1);
UIT0=subst(subst(UIT0,XX0,Y0),XX1,Y1);
UIT1=subst(subst(UIT1,XX0,Y0),XX1,Y1);
return(0)
}

k=3;
print("The number of compositions is: ",k)
for(i=1,k,compose())

een=1.0*(a*b)^(-1);
print("The inverse of ab is approximately: ",een)

/*
print("Here is UIT0:")
print(UIT0)
print("Here is UIT1:")
print(UIT1)
*/

HR(angle)=
{
	max(\
		abs(subst(subst(UIT0,X0,cos(angle)),X1,sin(angle))),\
		abs(subst(subst(UIT1,X0,cos(angle)),X1,sin(angle)))\
	)^(4^(-k))
}

twee=intnum(angle=0,2*Pi,HR(angle)^(-2));
print("The volume is approximately: ", twee)

print("a*b*Volume = ",a*b*twee)
