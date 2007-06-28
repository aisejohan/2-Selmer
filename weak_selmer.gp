\r n_linear_forms.gp.run

/* See equations in 1-selmer.gp */
has_points_at_p(a,b,d,u,v,p)=
{
local(m);
m=matrix(4,2);
m[1,1]=u*v*d;
m[1,2]=0;
m[2,1]=u^(-1)*d;
m[2,2]=u^(-1)*d*(-a);
m[3,1]=v^(-1)*d;
m[3,2]=v^(-1)*d*b;
m[4,1]=0;
m[4,2]=1;
return(n_linear_forms(m,p))
}

/* See equations in 1-selmer.gp */
has_points_at_2(a,b,d,u,v,k)=
{
local(m);
m=matrix(4,2);
m[1,1]=u*v*d;
m[1,2]=0;
m[2,1]=u^(-1)*d;
m[2,2]=u^(-1)*d*(-a);
m[3,1]=v^(-1)*d;
m[3,2]=v^(-1)*d*b;
m[4,1]=0;
m[4,2]=1;
return(stupid_solve_2(m,k)[1])
}

/* See equations in 1-selmer.gp */
has_points_at_infty(a,b,d,u,v)=
{
local(top,bottom);
bottom=0;
if(u>0,bottom=max(bottom,-b*u*v*d));
if(v>0,bottom=max(bottom,a*u*v*d));
if(u<0 && v>0,top=-b*u*v*d; return((top > bottom)));
if(u>0 && v<0,top=a*u*v*d; return((top > bottom)));
if(u<0 && v<0,top=min(-b*u*v*d,a*u*v*d); return((top > bottom)));
return(1)
}

/* Test if the ones having local points at p
form a linear subspace. */
test_linear_subspace(a,b,d,p)=
{
local(l,i,n,u,v,X,T);
if(p==2,error("Prime should not be 2."));
/* Always use the first nonsquare! */
n=2;
l=[];
while(issquare(Mod(n,p)),n++);
forvec(X=[[0,1],[0,1],[0,1],[0,1]],
	u=n^X[1]*p^X[2];
	v=n^X[3]*p^X[4];
	if(has_points_at_p(a,b,d,u,v,p),
		l=concat(l,[[Mod(X[1],2),Mod(X[2],2),Mod(X[3],2),Mod(X[4],2)]])
	)
);
T=matrix(4,matsize(l)[2],i,j,l[j][i]);
print("The matrix has rank ",matrank(T));
print("The matrix has size ",matsize(T));
print("The columns correspond to the ones that do have local points.");
return(T)
}

/* Test if the ones having local points at p
form a linear subspace. */
test_linear_subspace_2(a,b,d)=
{
local(g,l,ln,lr,i,ii,n,s,r,u,v,X,T,TT,k,skip);
if(valuation(a,2),
	g=2^valuation(a,2);
	a=a/g;
	b=b/g;
	d=d*g
);	
l=[];
ln=[];
lr=[];
r=0;
s=0;
n=0;
forvec(X=[[0,1],[0,1],[0,1],[0,1],[0,1],[0,1]],
	u=(-1)^X[1]*5^X[2]*2^X[3];
	v=(-1)^X[4]*5^X[5]*2^X[6];
	if(hilbert(v*(a+b)*b, u*a*(a+b), 2)==1 &&
	   hilbert(v, -u*d*a, 2)==1 &&
	   hilbert(u, v*d*b, 2)==1,
		if(has_points_at_2(a,b,d,u,v,10),
			s++;
			l=concat(l,
				[[Mod(X[1],2),Mod(X[2],2),Mod(X[3],2),
				Mod(X[4],2),Mod(X[5],2),Mod(X[6],2)]])
		,
			r++;
			lr=concat(lr,[[X[1],X[2],X[3],X[4],X[5],X[6]]])
		)
	,
		n++;
		ln=concat(ln,
				[[Mod(X[1],2),Mod(X[2],2),Mod(X[3],2),
				Mod(X[4],2),Mod(X[5],2),Mod(X[6],2)]])
	)
);
T=matrix(6,s,i,j,l[j][i]);
print("The matrix has rank ",matrank(T));
print("The number of vectors we found with local points is ",s);
print("The number of vectors provably having NO local points is ",n);
print("Hence the rank can be at most: ",floor(log(2^6-n)/log(2)));
k=12;
ii=1;
while(matrank(T)<3,
	X=lr[ii];
	TT=matrix(6,s+1);
	for(i=1,6,for(j=1,s,TT[i,j]=T[i,j]));
	for(i=1,6,TT[i,s+1]=Mod(X[i],2));
	if(matrank(T) < matrank(TT),
		u=(-1)^X[1]*5^X[2]*2^X[3];
		v=(-1)^X[4]*5^X[5]*2^X[6];
print(a," ",b," ",d," ",u," ",v," ",k);
		if(has_points_at_2(a,b,d,u,v,k),T=TT;s=s+1)
	);
	if(ii<r,ii++,ii=1;k=k+2)
);
print("The columns correspond to the ones that do have local points.");
return(T)
}

/* Test if the ones having local points at the reals
form a linear subspace. */
test_linear_infty(a,b,d)=
{
local(l,i,n,u,v,X,T);
/* Always use -1 as nonsquare! */
n=-1;
l=[];
forvec(X=[[0,1],[0,1]],
	u=n^X[1];
	v=n^X[2];
	if(has_points_at_infty(a,b,d,u,v),
		l=concat(l,[[Mod(X[1],2),Mod(X[2],2)]])
	)
);
T=matrix(2,matsize(l)[2],i,j,l[j][i]);
print("The matrix has rank ",matrank(T));
print("The matrix has size ",matsize(T));
print("The columns correspond to the ones that do have local points.");
return(T)
}

linear_conditions(a,b,d,p)=
{
local(l,i,j,n,u,v,X,T);
if(p==2,error("Prime should not be 2."));
/* Always use the first nonsquare! */
n=2;
l=[];
while(issquare(Mod(n,p)),n++);
forvec(X=[[0,1],[0,1],[0,1],[0,1]],
	u=n^X[1]*p^X[2];
	v=n^X[3]*p^X[4];
	if(has_points_at_p(a,b,d,u,v,p),
		l=concat(l,[[Mod(X[1],2),Mod(X[2],2),Mod(X[3],2),Mod(X[4],2)]])
	)
);
T=matrix(4,matsize(l)[2],i,j,l[j][i]);
if(!(matsize(l)[2]==matrank(T)^2),error("Not a linear subspace!"));
T=Mod(matker(T~),2)~;
return(T)
}

linear_conditions_2(a,b,d)=
{
local(g,l,ln,lr,i,ii,n,s,r,u,v,X,T,TT,k,skip);
if(valuation(a,2),
	g=2^valuation(a,2);
	a=a/g;
	b=b/g;
	d=d*g
);	
l=[];
ln=[];
lr=[];
r=0;
s=0;
n=0;
forvec(X=[[0,1],[0,1],[0,1],[0,1],[0,1],[0,1]],
	u=(-1)^X[1]*5^X[2]*2^X[3];
	v=(-1)^X[4]*5^X[5]*2^X[6];
	if(hilbert(v*(a+b)*b, u*a*(a+b), 2)==1 &&
	   hilbert(v, -u*d*a, 2)==1 &&
	   hilbert(u, v*d*b, 2)==1,
		if(has_points_at_2(a,b,d,u,v,10),
			s++;
			l=concat(l,
				[[Mod(X[1],2),Mod(X[2],2),Mod(X[3],2),
				Mod(X[4],2),Mod(X[5],2),Mod(X[6],2)]])
		,
			r++;
			lr=concat(lr,[[X[1],X[2],X[3],X[4],X[5],X[6]]])
		)
	,
		n++;
		ln=concat(ln,
				[[Mod(X[1],2),Mod(X[2],2),Mod(X[3],2),
				Mod(X[4],2),Mod(X[5],2),Mod(X[6],2)]])
	)
);
T=matrix(6,s,i,j,l[j][i]);
k=12;
ii=1;
while(matrank(T)<3,
	X=lr[ii];
	TT=matrix(6,s+1);
	for(i=1,6,for(j=1,s,TT[i,j]=T[i,j]));
	for(i=1,6,TT[i,s+1]=Mod(X[i],2));
	if(matrank(T) < matrank(TT),
		u=(-1)^X[1]*5^X[2]*2^X[3];
		v=(-1)^X[4]*5^X[5]*2^X[6];
print(a," ",b," ",d," ",u," ",v," ",k);
		if(has_points_at_2(a,b,d,u,v,k),T=TT;s=s+1)
	);
	if(ii<r,ii++,ii=1;k=k+2)
);
T=Mod(matker(T~),2)~;
return(T)
}

linear_conditions_infty(a,b,d)=
{
local(l,i,n,u,v,X,T);
/* Always use -1 as nonsquare! */
n=-1;
l=[];
forvec(X=[[0,1],[0,1]],
	u=n^X[1];
	v=n^X[2];
	if(has_points_at_infty(a,b,d,u,v),
		l=concat(l,[[Mod(X[1],2),Mod(X[2],2)]])
	)
);
T=matrix(2,matsize(l)[2],i,j,l[j][i]);
T=Mod(matker(T~),2)~;
return(T)
}

global_linear_conditions(a,b,d)=
{
local(s,l,N,T,UIT,i,j,e,p,t,tmp,powern,powerp,power5);
N=-abs(a*b*(a+b)*d);
l=factor(N);
s=matsize(l)[1];
UIT=matrix(2*s,2*s);
/* Real numbers. */
T=linear_conditions_infty(a,b,d);
for(j=1,s,
	if(sign(l[j,1])==1,
		powern=Mod(0,2)
	,
		powern=Mod(1,2)
	);
	t=[powern,Mod(0,2)]~;
	tmp=T*t;
	UIT[j,1]=tmp[1];
	t=[Mod(0,2),powern]~;
	tmp=T*t;
	UIT[s+j,1]=tmp[1]
);
/* 2-adic */
T=linear_conditions_2(a,b,d);
for(j=1,s,
	e=valuation(l[j,1],2);
	powerp=e;
	powern=1-issquare(Mod(2^(-e)*l[j,1],4));
	power5=1-issquare(Mod(2^(-e)*(-1)^(powern)*l[j,1],8));
if(!issquare(Mod(l[j,1]*(-1)^powern*5^power5*2^powerp,256)),error("FUCK"));
	powerp=Mod(powerp,2);
	powern=Mod(powern,2);
	power5=Mod(power5,2);
	t=[powern,power5,powerp,Mod(0,2),Mod(0,2),Mod(0,2)]~;
	tmp=T*t;
	UIT[j,2]=tmp[1];
	UIT[j,3]=tmp[2];
	UIT[j,4]=tmp[3];
	t=[Mod(0,2),Mod(0,2),Mod(0,2),powern,power5,powerp]~;
	tmp=T*t;
	UIT[s+j,2]=tmp[1];
	UIT[s+j,3]=tmp[2];
	UIT[s+j,4]=tmp[3]
);
/* Odd primes. */
for(i=3,s,
	p=l[i,1];
	T=linear_conditions(a,b,d,p);
	for(j=1,s,
		e=valuation(l[j,1],p);
		powerp=Mod(e,2);
		powern=Mod(1-issquare(Mod(p^(-e)*l[j,1],p)),2);
		t=[powern,powerp,Mod(0,2),Mod(0,2)]~;
		tmp=T*t;
		UIT[j,2*i-1]=tmp[1];
		UIT[j,2*i]=tmp[2];
		t=[Mod(0,2),Mod(0,2),powern,powerp]~;
		tmp=T*t;
		UIT[s+j,2*i-1]=tmp[1];
		UIT[s+j,2*i]=tmp[2];
	);
);
return(matsize(UIT)[1]-matrank(UIT))
}
