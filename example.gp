find_l(a,b,d)=
{
local(N,l);
N=-abs(a*b*(a+b)*d);
l=factor(N)[,1]~;
return(l)
}

find_triv(a,b,d)=
{
local(triv1,triv2,triv3,l,s,i,p);
l=find_l(a,b,d);
s=matsize(l)[2];
triv1=vector(2*s);
triv2=vector(2*s);
triv3=vector(2*s);
/* infinity */
if(-a*d < 0, triv1[1]=1);
if(a*(a+b) < 0, triv2[1]=1);
if(-(a+b)*d < 0, triv3[1]=1);
if(b*d < 0, triv1[s+1]=1);
if((a+b)*d < 0, triv2[s+1]=1);
if((a+b)*b < 0, triv3[s+1]=1);
/* powers of primes */
for(i=2,s,
	p=l[i];
	if(Mod(valuation(-a*d,p),2), triv1[i]=1);
	if(Mod(valuation(a*(a+b),p),2), triv2[i]=1);
	if(Mod(valuation(-(a+b)*d,p),2), triv3[i]=1);
	if(Mod(valuation(b*d,p),2), triv1[s+i]=1);
	if(Mod(valuation((a+b)*d,p),2), triv2[s+i]=1);
	if(Mod(valuation((a+b)*b,p),2), triv3[s+i]=1);
);
return([triv1,triv2,triv3])
}

uv_vector(a,b,d,vv)=
{
local(i,u,v,l);
l=find_l(a,b,d);
s=matsize(l)[2];
if(matsize(vv)[2] == 2*s,,
	error("Wrong sizes.")
);
u=1;
v=1;
for(i=1,s,
	u=u*l[i]^vv[i];
	v=v*l[i]^vv[s+i]
);
return([u,v])
}

quadrics(a,b,d,uv)=
{
local(u,v,G1,G2,tmp);
u=uv[1];
v=uv[2];
tmp=[0,u,-v,(a+b)*d];
G1=matdiagonal(tmp);
tmp=[-(a+b)*(u*v)^(-1),b*u,a*v,0];
G2=matdiagonal(tmp);
return([G1,G2])
}

/* The point will be on the curve
	d*y^2 = x*(x-a)*(x+b)	*/
associated_point(a,b,d,uv,oplossing)=
{
local(x,y,u,v);
u = uv[1];
v = uv[2];
x = oplossing[1]^2/(u*v*d*oplossing[4]^2);
y = oplossing[1]*oplossing[2]*oplossing[3]/(oplossing[4]^3*d^2);
return([x,y])
}


try_it(aa,bb,dd,uv,GG,B)=
{
local(i,j);
GG[1]=GG[1]/content(GG[1]);
GG[2]=GG[2]/content(GG[2]);
for(i=-B,B,
	for(j=1,B,
		if(gcd(i,j) == 1,,
			www=Qfsolve(i*GG[1]+j*GG[2])~;
			if((www*GG[1])*(www~),,
				print (
				"Here (i,j)=(",i,",",j,")",
				"\tand www=",www,
				"\tand point=",
				associated_point(aa,bb,dd,uv,www));
				return(0);
			)
		)
	)
);
return(1)
}

disc_at_p(G,P)=
{
local(M,lattice,Q,v1,v2);
P = P/content(P);
M=matrix(1,4,i,j,(P*G)[j]);
lattice=matkerint(M);
M=matrix(4,3);
M[,1]=P~;
M[,2]=lattice[,1];
M[,3]=lattice[,2];
if(matrank(M) == 3,
	v1=lattice[,1];
	v2=lattice[,2]
,
	error("FIXME")
);
Q=matrix(2,2);
Q[1,1]=(v1~)*G*v1;
Q[1,2]=(v1~)*G*v2;
Q[2,1]=(v2~)*G*v1;
Q[2,2]=(v2~)*G*v2;
return(Qfsolve(Q))
}

maak_lijst(B)=
{
local(lijst,i1,i2,i3,i4);
lijst=[];
i1=0;
while(i1 <= B,
	if(i1,
		i2=-B
	,
		i2=0
	);
	while( i2 <= B,
		if(i1 || i2,
			i3=-B
		,
			i3=0
		);
		while(i3 <= B,
			if(i1 || i2 || i3,
				i4=-B
			,
				i4=1
			);
			while(i4 <= B,
				if(gcd([i1,i2,i3,i4]) == 1,
					lijst=concat(lijst,[[i1,i2,i3,i4]])
				);
				i4++
			);
			i3++
		);
		i2++
	);
	i1++
);
return(lijst)
}

find_smaller(G,P,lijst)=
{
local(t,H,X);
P=P/content(P);
H=sum(i=1,4,P[i]^2);
for(i=1,matsize(lijst)[2],
	X=lijst[i];
	t=-2*X*G*P~/(X*G*X~);
	X=P+t*X;
	X=X/content(X);
	print(X,"    ",X*G*X~);
	if(H > sum(j=1,4,X[j]^2),return(X));
);
return(0)
}

aa=1;
bb=2;
dd=115;

UU=global_linear_conditions(aa,bb,dd);
KERN=lift(matker(UU~))~;
rank=matsize(KERN)[1];
print("The rank of the 2-selmer group is: ",rank);
l=find_l(aa,bb,dd);
trivials=find_triv(aa,bb,dd);
tmp=vector(rank,i,[0,1]);
forvec(X=tmp,\
{
	test=lift(Mod(X*KERN,2));
	if( test == 0 || 
	test == trivials[1] || 
	test == trivials[2] ||
	test == trivials[3],,
		uv=uv_vector(aa,bb,dd,test);
		GG=quadrics(aa,bb,dd,uv);
		if(try_it(aa,bb,dd,uv,GG,50),
			print("Does not work when X=",X)
		)
	)
})
