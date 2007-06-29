n_helper(t,v)=
{
	if(t,
		if(t[1]*v[1]+t[2]*v[2],
			return([0,0])
		,
			return(t)
		)
	,
		if(v[2],
			return([-v[2],v[1]])
		,
			return([0,0])
		)
	)
}

n_linear_forms_mod_p(m,p)=
{
local(i,j,e,n,t,l,h,s);
n=matsize(m)[1];
t=vector(2);
l=[];
s=1;
for(i=1,n,
	e=valuation(m[i,],p);
	m[i,]=Mod(p^(-e)*m[i,],p);
	if(Mod(e,2),
		t=n_helper(t,m[i,]);
		if(!t,
			return([0])
		)
	,
		l=concat(l,[i]);
	)
);
n=matsize(l)[2];
for(i=1,n,
	for(j=i+1,n,
		if(!(m[l[i],1]*m[l[j],2]-m[l[i],2]*m[l[j],1]),
			if(
			(
			m[l[i],1] &&
			!issquare(m[l[j],1]/m[l[i],1])
			) ||
			(
			m[l[i],2] &&
			!issquare(m[l[j],2]/m[l[i],2])
			),
				t=n_helper(t,m[l[i],]);
				if(!t,
					return([0])
				)
			)
		)
	)
);
if(t,
	for(i=1,n,
		h=t[1]*m[l[i],1]+t[2]*m[l[i],2];
		if(h,
			if(!issquare(h),
				t=t*h
			);
			s=2;
			i=n
		)
	);
	for(i=1,n,
		if(!issquare(t[1]*m[l[i],1]+t[2]*m[l[i],2]),
				return([0])
		)
	)
);
return([s,t])
}

n_return(t,m,p)=
{
local(b);
b=[t[1],0;t[2],p];
return(m*b)
}

/* If s=1 looks for a solution not of the form [0,*] mod p.*/
/* If s=2 looks for a solution of the form [1,*] mod p. */
n_linear_forms_recursive(m,p,s)=
{
local(t,u,e,b,i,j,n,c,ss);
if(!(s==1) && !(s==2),
	error("Wrong s!")
);
n=matsize(m)[1];
i=0;
while(i<n,
	i++;
	if(m[i,1],
		e=valuation(m[i,1],p);
		if(Mod(e,2),
			i=n+1
		,
			if(!issquare(Mod(p^(-e)*m[i,1],p)),
				i=n+1
			)
		)
	)
);
if(i==n,
	return(1)
);
u=n_linear_forms_mod_p(m,p);
if(u[1],
	if(u[2],
		if(s==2 && u[2]==2 && !issquare(u[2][1]),
			return(0)
		);
		if(u[2]==2,
			s=2
		);
		t=[lift(u[2][1]),lift(u[2][2])];
		return(n_linear_forms_recursive(n_return(t,m,p),p,s))
	,
		/* In this case all the exponents are even! */
		for(i=1,n,
			e=valuation(m[i,],p);
			if(Mod(e,2),
				error("Quel horreur!")
			);
			m[i,]=p^(-e)*m[i,]
		);
		for(i=0,p-1,
			t=[1,i];
			j=1;
			ss=s;
			while(j<n+1 && !c=Mod(m[j,]*t~,p),
				j++
			);
			if(!issquare(c),
				if(s==2,
					j=n+2
				,
					t=t*lift(c);
					ss=2
				)
			);
			while(j<n+1 && issquare(Mod(m[j,]*t~,p)),
				j++
			);
			if(j==n+1,
				if(n_linear_forms_recursive(n_return(t,m,p),p,ss),
					return(1)
				)
			)
		);
		return(0)
	)
,
	return(0)
)
}

n_linear_forms(m,p)=
{
local(i,j,e,n,t,c);
/* Look for trivial solutions first. */
n=matsize(m)[1];
for(i=1,n,\
	if(!m[i,],
		error("Equation is zero. Not allowed!")
	);
	t=[-m[i,2],m[i,1]];
	j=1;
	while(j<n+1 && !c=m[j,]*t~,
		j++
	);
	t=t*c;
	while(j<n+1,
		if(c=m[j,]*t~,
			e=valuation(c,p);
			if(Mod(e,2),
				j=n+2
			,
				if(!issquare(Mod(p^(-e)*c,p)),
					j=n+2
				)
			)
		);
		j++
	);
	if(j==n+1,return(1))
);
if(n_linear_forms_recursive(m,p,1),
	return(1)
,
	return(n_linear_forms_recursive(m*[0,p;1,0],p,1))
)
}

stupid_solve(m,p,k)=
{
local(X,Y,t,i,j,e,n,c);
n=matsize(m)[1];
for(i=0,k,
	X=p^i;
	for(Y=0,p^(k-i)-1,
		t=[X,Y];
		if(valuation([X,Y],p)==0,
			j=1;
			while(j<n+1 && !c=m[j,]*t~,
				j++
			);
			t=t*c;
			while(j<n+1,
				if(c=m[j,]*t~,
					e=valuation(c,p);
					if(Mod(e,2),
						j=n+2
					,
						if(!issquare(Mod(p^(-e)*c,p)),
							j=n+2
						)
					)
				);
				j++
			);
			if(j==n+1,
				return([1,t])
			)
		)
	)
);
return([0,"No solution"])
}

test_function(a,b,p,k)=
{
local(l,m,mm,n,u,v,d,X,smart,stupid);
if(p==2,error("Prime should not be 2."));
n=2;
while(issquare(Mod(n,p)),n++);
l=vector(4);
l[1]=1;
l[2]=n;
l[3]=p;
l[4]=n*p;
m=matrix(4,2);
forvec(X=[[1,4],[1,4],[1,4]],
	u=l[X[1]];
	v=l[X[2]];
	d=l[X[3]];
	m[1,1]=u*v*d;
	m[1,2]=0;
	m[2,1]=u^(-1)*d;
	m[2,2]=u^(-1)*d*(-a);
	m[3,1]=v^(-1)*d;
	m[3,2]=v^(-1)*d*b;
	m[4,1]=0;
	m[4,2]=1;
	smart=n_linear_forms(m,p);
	stupid=stupid_solve(m,p,k);
	if(!(smart==stupid[1]),
		print("AAAAAAAAAAAARRRRRRRRRRRRRRGGGGGGGGGGGGGHHHHHHHHH");
		print(X);
		print(m);
		print(stupid[2])
	)
)
}

test_function_2(p)=
{
local(l,m,mm,n,u,v,d,X,smart,stupid);
if(p==2,
	error("Prime should not be 2.")
);
n=2;
while(issquare(Mod(n,p)),
	n++
);
l=vector(4);
l[1]=1;
l[2]=n;
l[3]=p;
l[4]=n*p;
m=matrix(3,2);
forvec(X=[[1,4],[1,4]],
	u=l[X[1]];
	v=l[X[2]];
	m[1,1]=u;
	m[1,2]=0;
	m[2,1]=u*v;
	m[2,2]=u*v;
	m[3,1]=0;
	m[3,2]=v;
	smart=n_linear_forms(m,p);
	stupid=1;
	if(X[1]==2 && X[2] == 3,
		stupid=0
	);
	if(X[1]==2 && X[2] == 4,
		stupid=0
	);
	if(X[1]==3 && X[2] == 2,
		stupid=0
	);
	if(X[1]==4 && X[2] == 2,
		stupid=0
	);
	if(issquare(Mod(-1,p)),
		if(X[1]==3 && X[2] ==4,
			stupid=0
		);
		if(X[1]==4 && X[2] ==3,
			stupid=0
		)
	,
		if(X[1]==3 && X[2] ==3,
			stupid=0
		);
		if(X[1]==4 && X[2] ==4,
			stupid=0
		)
	);
	if(!(smart==stupid),
		print("AAAAAAAAAAAARRRRRRRRRRRRRRGGGGGGGGGGGGGHHHHHHHHH");
		print(X);
		print(m);
	);
	stupid=hilbert(u,v,p);
	if(stupid==-1,
		stupid=0
	);
	if(!(smart==stupid),
		print("AAAAAAAAAAAARRRRRRRRRRRRRRGGGGGGGGGGGGGHHHHHHHHH");
		print(X);
		print(m);
	);
)
}

stupid_solve_2(m,k)=
{
local(X,Y,t,i,j,e,n,c);
n=matsize(m)[1];
/* Look for trivial solutions; and scale equations. */
for(i=1,n,\
	if(!m[i,],
		error("Equation is zero. Not allowed!")
	);
	e=valuation(m[i,],2);
	if(Mod(e,2),
		m[i,]=2^(-e+1)*m[i,]
	,
		m[i,]=2^(-e)*m[i,]
	);
	t=[-m[i,2],m[i,1]];
	j=1;
	while(j<n+1 && !c=m[j,]*t~,
		j++
	);
	t=t*c;
	while(j<n+1,
		if(c=m[j,]*t~,
			e=valuation(c,2);
			if(Mod(e,2),
				j=n+2
			,
				if(!issquare(Mod(2^(-e)*c,16)),
					j=n+2
				)
			)
		);
		j++
	);
	if(j==n+1,
		return([1,t])
	)
);
for(i=0,k,
	X=2^i;
	for(Y=0,2^(k-i)-1,
		t=[X,Y];
		if(valuation([X,Y],2)==0,
			j=1;
			while(j<n+1 && !c=m[j,]*t~,
				j++
			);
			t=t*c;
			while(j<n+1,
				if(c=m[j,]*t~,
					e=valuation(c,2);
					if(Mod(e,2),
						j=n+2
					,
						if(!issquare(Mod(2^(-e)*c,16)),
							j=n+2
						)
					)
				);
				j++
			);
			if(j==n+1,return([1,t]))
		)
	)
);
return([0,"No solution"])
}

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
if(u>0,
	bottom=max(bottom,-b*u*v*d)
);
if(v>0,
	bottom=max(bottom,a*u*v*d)
);
if(u<0 && v>0,
	top=-b*u*v*d;
	return((top > bottom))
);
if(u>0 && v<0,
	top=a*u*v*d;
	return((top > bottom))
);
if(u<0 && v<0,
	top=min(-b*u*v*d,a*u*v*d);
	return((top > bottom))
);
return(1)
}

/* Test if the ones having local points at p
form a linear subspace. */
test_linear_subspace(a,b,d,p)=
{
local(l,i,n,u,v,X,T);
if(p==2,
	error("Prime should not be 2.")
);
/* Always use the first nonsquare! */
n=2;
l=[];
while(issquare(Mod(n,p)),
	n++
);
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
	for(i=1,6,
		for(j=1,s,
			TT[i,j]=T[i,j]
		)
	);
	for(i=1,6,
		TT[i,s+1]=Mod(X[i],2)
	);
	if(matrank(T) < matrank(TT),
		u=(-1)^X[1]*5^X[2]*2^X[3];
		v=(-1)^X[4]*5^X[5]*2^X[6];
		print(a," ",b," ",d," ",u," ",v," ",k);
		if(has_points_at_2(a,b,d,u,v,k),
			T=TT;
			s=s+1
		)
	);
	if(ii<r,
		ii++
	,
		ii=1;
		k=k+2
	)
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
while(issquare(Mod(n,p)),
	n++
);
forvec(X=[[0,1],[0,1],[0,1],[0,1]],
	u=n^X[1]*p^X[2];
	v=n^X[3]*p^X[4];
	if(has_points_at_p(a,b,d,u,v,p),
		l=concat(l,[[Mod(X[1],2),Mod(X[2],2),Mod(X[3],2),Mod(X[4],2)]])
	)
);
T=matrix(4,matsize(l)[2],i,j,l[j][i]);
if(!(matsize(l)[2]==matrank(T)^2),
	error("Not a linear subspace!")
);
T=Mod(matker(T~),2)~;
return(T)
}

set_of_linear_conditions(a,b,p)=
{
local(n,uit);
n=2;
while(issquare(Mod(n,p)),
	n++
);
uit=vector(4);
uit[1]=linear_conditions(a,b,1,p);
uit[2]=linear_conditions(a,b,n,p);
uit[3]=linear_conditions(a,b,p,p);
uit[4]=linear_conditions(a,b,p*n,p);
return(uit)
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
	for(i=1,6,
		for(j=1,s,
			TT[i,j]=T[i,j]
		)
	);
	for(i=1,6,
		TT[i,s+1]=Mod(X[i],2)
	);
	if(matrank(T) < matrank(TT),
		u=(-1)^X[1]*5^X[2]*2^X[3];
		v=(-1)^X[4]*5^X[5]*2^X[6];
		print(a," ",b," ",d," ",u," ",v," ",k);
		if(has_points_at_2(a,b,d,u,v,k),
			T=TT;
			s=s+1
		)
	);
	if(ii<r,ii++,
		ii=1;
		k=k+2
	)
);
T=Mod(matker(T~),2)~;
return(T)
}

set_of_linear_conditions_2(a,b)=
{
local(n,uit);
uit=vector(8);
uit[1]=linear_conditions_2(a,b,1);
uit[2]=linear_conditions_2(a,b,-1);
uit[3]=linear_conditions_2(a,b,5);
uit[4]=linear_conditions_2(a,b,-5);
uit[5]=linear_conditions_2(a,b,2);
uit[6]=linear_conditions_2(a,b,-2);
uit[7]=linear_conditions_2(a,b,10);
uit[8]=linear_conditions_2(a,b,-10);
return(uit)
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

set_of_linear_conditions_infty(a,b)=
{
local(n,uit);
uit=vector(2);
uit[1]=linear_conditions_infty(a,b,1);
uit[2]=linear_conditions_infty(a,b,-1);
return(uit)
}

/* Create the list of sets of local linear conditions. */
j_init(a,b,j_lim)=
{
local(lijst,i);
if(j_lim<3,
	error("Come on!")
);
lijst=vector(j_lim+1);
lijst[1]=set_of_linear_conditions_infty(a,b);
lijst[2]=set_of_linear_conditions_2(a,b);
for(i=3,j_lim+1,
	lijst[i]=set_of_linear_conditions(a,b,prime(i-1))
);
return(lijst)
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
	if(!issquare(Mod(l[j,1]*(-1)^powern*5^power5*2^powerp,256)),
		error("Quel horreur!")
	);
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
return(UIT);
}

global_linear_conditions_lijst(a,b,d,lijst)=
{
local(s,l,N,T,UIT,i,j,e,p,t,tmp,powern,powerp,power5);
N=-abs(a*b*(a+b)*d);
l=factor(N);
s=matsize(l)[1];
UIT=matrix(2*s,2*s);
/* Real numbers. */
if(sign(d) == 1,
	T=lijst[1][1]
,
	T=lijst[1][2]
);
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
powerp=valuation(d,2);
powern=1-issquare(Mod(2^(-powerp)*d,4));
power5=1-issquare(Mod(2^(-powerp)*(-1)^powern*d,8));
powerp=lift(Mod(powerp,2));
T=lijst[2][1+powern+2*power5+4*powerp];
for(j=1,s,
	e=valuation(l[j,1],2);
	powerp=e;
	powern=1-issquare(Mod(2^(-e)*l[j,1],4));
	power5=1-issquare(Mod(2^(-e)*(-1)^(powern)*l[j,1],8));
	if(!issquare(Mod(l[j,1]*(-1)^powern*5^power5*2^powerp,256)),
		error("Quel horreur!")
	);
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
	j=primepi(p);
	if(p == prime(j),,
		error("Quel horreur!")
	);
	powerp=valuation(d,p);
	powern=1-issquare(Mod(p^(-powerp)*d,p));
	powerp=lift(Mod(powerp,2));
	T=lijst[j+1][1+powern+2*powerp];
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
return(UIT);
}

test_global_lin_con_twice(aa,bb)=
{
local(d,u1,u2,maxpi,lijst);
maxpi=50;
maxpi=max(maxpi,primepi(abs(aa)));
maxpi=max(maxpi,primepi(abs(bb)));
maxpi=max(maxpi,primepi(abs(aa+bb)));
lijst=j_init(aa,bb,maxpi+1);
print(maxpi+1);
d=-prime(45);
while(d <= prime(45),
	print(d);
	if(d==0,
		d++
	);
	u1=global_linear_conditions(aa,bb,d);
	u2=global_linear_conditions_lijst(aa,bb,d,lijst);
	if(u1==u2,
		print("YES YES YES YES")
	,
		print("NO NO NO NO")
	);
	d++
);
return(0)
}

rank_2_selmer_group(a,b,d)=
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
	if(!issquare(Mod(l[j,1]*(-1)^powern*5^power5*2^powerp,256)),
		error("Quel horreur!")
	);
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

rank_2_selmer_group_lijst(a,b,d,lijst)=
{
local(s,l,N,T,UIT,i,j,e,p,t,tmp,powern,powerp,power5);
N=-abs(a*b*(a+b)*d);
l=factor(N);
s=matsize(l)[1];
UIT=matrix(2*s,2*s);
/* Real numbers. */
if(sign(d) == 1,
	T=lijst[1][1]
,
	T=lijst[1][2]
);
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
powerp=valuation(d,2);
powern=1-issquare(Mod(2^(-powerp)*d,4));
power5=1-issquare(Mod(2^(-powerp)*(-1)^powern*d,8));
powerp=lift(Mod(powerp,2));
T=lijst[2][1+powern+2*power5+4*powerp];
for(j=1,s,
	e=valuation(l[j,1],2);
	powerp=e;
	powern=1-issquare(Mod(2^(-e)*l[j,1],4));
	power5=1-issquare(Mod(2^(-e)*(-1)^(powern)*l[j,1],8));
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
	j=primepi(p);
	powerp=valuation(d,p);
	powern=1-issquare(Mod(p^(-powerp)*d,p));
	powerp=lift(Mod(powerp,2));
	T=lijst[j+1][1+powern+2*powerp];
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

j_compute(j_a,j_b,j_d)=
{
	local(j_return);
	j_return=rank_2_selmer_group(j_a,j_b,j_d);
	if(j_return > 10,
		print("Here the 2-Selmer rank is: ",j_return);
		write("BIGSELMER",
		"a=",j_a,"; b=",j_b,"; d="j_d,"; selmer=",j_return,";")
	);
	return(j_return)
}

j_compute_lijst(j_a,j_b,j_d,lijst)=
{
	local(j_return);
	j_return=rank_2_selmer_group_lijst(j_a,j_b,j_d,lijst);
	if(j_return > 10,
		print("Here the 2-Selmer rank is: ",j_return);
		write("BIGSELMER",
		"a=",j_a,"; b=",j_b,"; d="j_d,"; selmer=",j_return,";")
	);
	return(j_return)
}

j_setup(j_a,j_b,j_bound,j_lim)=
{
	local(l,len,j_prime,j_rank);
	l=[];
	len=0;
	j_rank=j_compute(j_a,j_b,-1);
	if(j_rank >= j_bound,
		len++;
		l=concat(l,[-1])
	);
	j_prime=2;
	while(j_prime < j_lim,
		j_rank=j_compute(j_a,j_b,j_prime);
		if(j_rank >= j_bound,
			len++;
			l=concat(l,[j_prime])
		);
		j_prime=nextprime(j_prime+1)
	);
	return(l)
}

j_setup_lijst(j_a,j_b,j_bound,j_lim,lijst)=
{
	local(l,len,j_index,j_prime,j_rank);
	l=[];
	len=0;
	j_rank=j_compute_lijst(j_a,j_b,-1,lijst);
	if(j_rank >= j_bound,
		len++;
		l=concat(l,[-1])
	);
	j_index=1;
	while(j_index <= j_lim,
		j_prime=prime(j_index);
		j_rank=j_compute_lijst(j_a,j_b,j_prime,lijst);
		if(j_rank >= j_bound,
			len++;
			l=concat(l,[j_prime])
		);
		j_index++
	);
	return(l)
}

/* Look for index smallest element q in the list such that target <= q.
   Returns length+1 if there is no such index. */
find_spot(j_lijst,j_target)=
{
	local(j_len,j_above,j_below,j_ii);
	j_len=matsize(j_lijst)[2];
	if(j_lijst[j_len] < j_target,
		return(j_len+1)
	,
		j_above=j_len
	);
	if(j_lijst[1] >= j_target,
		return(1)
	,
		j_below=1
	);
	while(j_below < j_above-1,
		j_ii=round((j_above+j_below)/2);
		if(j_lijst[j_ii] < j_target,
			j_below=j_ii
		,
			j_above=j_ii
		)
	);
	return(j_above)
}

/*
	Abbreviation: prev=previous
	l_primes is an ordered list of primes with possibly also -1 in the
	first spot.
	l_prev is an ordered list of square free integers each a product of
	exactly	the same number of factors where -1 counts as a factor.
*/
j_do_it_once(j_a,j_b,j_bound,l_primes,l_prev)=
{
	local(l,len,len_prev,len_primes,j_tmp,j_factors,j_j,j_prime,
		j_ok,j_n,j_spot,j_rank,j_got);
	l=[];
	len=0;
	len_primes=matsize(l_primes)[2];
	len_prev=matsize(l_prev)[2];
	if(len_prev == 0 || len_primes == 0,
		return(l)
	,
		j_got=matsize(factor(l_prev[1]))[1]
	);
	for(j_i=1,len_prev,
		print("Doing ",j_i," out of ",len_prev);
		j_tmp=l_prev[j_i];
		j_factors=factor(j_tmp)[,1];
		j_j=find_spot(l_primes,j_factors[j_got]+1);
		while(j_j <= len_primes,
			j_prime=l_primes[j_j];
			j_ok=1;
			for(j_k=1,j_got,
				j_n=j_tmp*j_prime/j_factors[j_k];
				j_spot=find_spot(l_prev,j_n);
				if(j_spot > len_prev || j_n != l_prev[j_spot],
					j_ok=0;
					j_k=j_got;
				)
			);
			if(j_ok,
				j_rank=j_compute(j_a,j_b,j_tmp*j_prime);
				if(j_rank >= j_bound,
					len++;
					l=concat(l,[j_tmp*j_prime])
				)
			);
			j_j++
		)
	);
	l=vecsort(l);
	return(l)
}

j_do_it_once_lijst(j_a,j_b,j_bound,l_primes,l_prev,lijst)=
{
	local(l,len,len_prev,len_primes,j_tmp,j_factors,j_j,j_prime,
		j_ok,j_n,j_spot,j_rank,j_got);
	l=[];
	len=0;
	len_primes=matsize(l_primes)[2];
	len_prev=matsize(l_prev)[2];
	if(len_prev == 0 || len_primes == 0,
		return(l)
	,
		j_got=matsize(factor(l_prev[1]))[1]
	);
	for(j_i=1,len_prev,
		print("Doing ",j_i," out of ",len_prev);
		j_tmp=l_prev[j_i];
		j_factors=factor(j_tmp)[,1];
		j_j=find_spot(l_primes,j_factors[j_got]+1);
		while(j_j <= len_primes,
			j_prime=l_primes[j_j];
			j_ok=1;
			for(j_k=1,j_got,
				j_n=j_tmp*j_prime/j_factors[j_k];
				j_spot=find_spot(l_prev,j_n);
				if(j_spot > len_prev || j_n != l_prev[j_spot],
					j_ok=0;
					j_k=j_got;
				)
			);
			if(j_ok,
				j_rank=j_compute_lijst(j_a,j_b,j_tmp*j_prime,lijst);
				if(j_rank >= j_bound,
					len++;
					l=concat(l,[j_tmp*j_prime])
				)
			);
			j_j++
		)
	);
	l=vecsort(l);
	return(l)
}

help()=
{
print("
n_helper(t,v)
n_linear_forms_mod_p(m,p)
n_return(t,m,p)
n_linear_forms_recursive(m,p,s)
n_linear_forms(m,p)
stupid_solve(m,p,k)
test_function(a,b,p,k)
test_function_2(p)
stupid_solve_2(m,k)
has_points_at_p(a,b,d,u,v,p)
has_points_at_2(a,b,d,u,v,k)
has_points_at_infty(a,b,d,u,v)
test_linear_subspace(a,b,d,p)
test_linear_subspace_2(a,b,d)
test_linear_infty(a,b,d)
linear_conditions(a,b,d,p)
set_of_linear_conditions(a,b,p)
linear_conditions_2(a,b,d)
set_of_linear_conditions_2(a,b)
linear_conditions_infty(a,b,d)
set_of_linear_conditions_infty(a,b)
j_init(a,b,j_lim)
global_linear_conditions(a,b,d)
global_linear_conditions_lijst(a,b,d,lijst)
test_global_lin_con_twice(aa,bb)
rank_2_selmer_group(a,b,d)
rank_2_selmer_group_lijst(a,b,d,lijst)
j_compute(j_a,j_b,j_d)
j_compute_lijst(j_a,j_b,j_d,lijst)
j_setup(j_a,j_b,j_bound,j_lim)
j_setup_lijst(j_a,j_b,j_bound,j_lim,lijst)
find_spot(j_lijst,j_target)
j_do_it_once(j_a,j_b,j_bound,l_primes,l_prev)
j_do_it_once_lijst(j_a,j_b,j_bound,l_primes,l_prev,lijst)
help()");
return(0)
}

