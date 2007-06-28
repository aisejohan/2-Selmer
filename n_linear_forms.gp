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
		if(!t,return([0]))
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
				if(!t,return([0]))
			)
		)
	)
);
if(t,
	for(i=1,n,
		h=t[1]*m[l[i],1]+t[2]*m[l[i],2];
		if(h,
			if(!issquare(h),t=t*h);
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
if(!(s==1) && !(s==2), error("Wrong s!"));
n=matsize(m)[1];
i=0;
while(i<n,
	i++;
	if(m[i,1],
		e=valuation(m[i,1],p);
		if(Mod(e,2),
			i=n+1
		,
			if(!issquare(Mod(p^(-e)*m[i,1],p)),i=n+1)
		)
	)
);
if(i==n,return(1));
u=n_linear_forms_mod_p(m,p);
if(u[1],
	if(u[2],
		if(s==2 && u[2]==2 && !issquare(u[2][1]),return(0));
		if(u[2]==2,s=2);
		t=[lift(u[2][1]),lift(u[2][2])];
		return(n_linear_forms_recursive(n_return(t,m,p),p,s))
	,
/* In this case all the exponents are even! */
		for(i=1,n,
			e=valuation(m[i,],p);
			if(Mod(e,2),error("Quel horreur!"));
			m[i,]=p^(-e)*m[i,]
		);
		for(i=0,p-1,
			t=[1,i];
			j=1;
			ss=s;
			while(j<n+1 && !c=Mod(m[j,]*t~,p),j++);
			if(!issquare(c),
				if(s==2,
					j=n+2
				,
					t=t*lift(c);
					ss=2
				)
			);
			while(j<n+1 && issquare(Mod(m[j,]*t~,p)),j++);
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
	if(!m[i,],error("Equation is zero. Not allowed!"));
	t=[-m[i,2],m[i,1]];
	j=1;
	while(j<n+1 && !c=m[j,]*t~,j++);
	t=t*c;
	while(j<n+1,
		if(c=m[j,]*t~,
			e=valuation(c,p);
			if(Mod(e,2),
				j=n+2
			,
				if(!issquare(Mod(p^(-e)*c,p)),j=n+2)
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
			while(j<n+1 && !c=m[j,]*t~,j++);
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
			if(j==n+1,return([1,t]))
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
if(p==2,error("Prime should not be 2."));
n=2;
while(issquare(Mod(n,p)),n++);
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
	if(X[1]==2 && X[2] == 3, stupid=0);
	if(X[1]==2 && X[2] == 4, stupid=0);
	if(X[1]==3 && X[2] == 2, stupid=0);
	if(X[1]==4 && X[2] == 2, stupid=0);
	if(issquare(Mod(-1,p)),
		if(X[1]==3 && X[2] ==4, stupid=0);
		if(X[1]==4 && X[2] ==3, stupid=0)
	,
		if(X[1]==3 && X[2] ==3, stupid=0);
		if(X[1]==4 && X[2] ==4, stupid=0)
	);
	if(!(smart==stupid),
		print("AAAAAAAAAAAARRRRRRRRRRRRRRGGGGGGGGGGGGGHHHHHHHHH");
		print(X);
		print(m);
	);
	stupid=hilbert(u,v,p);
	if(stupid==-1,stupid=0);
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
	if(!m[i,],error("Equation is zero. Not allowed!"));
	e=valuation(m[i,],2);
	if(Mod(e,2),
		m[i,]=2^(-e+1)*m[i,]
	,
		m[i,]=2^(-e)*m[i,]
	);
	t=[-m[i,2],m[i,1]];
	j=1;
	while(j<n+1 && !c=m[j,]*t~,j++);
	t=t*c;
	while(j<n+1,
		if(c=m[j,]*t~,
			e=valuation(c,2);
			if(Mod(e,2),
				j=n+2
			,
				if(!issquare(Mod(2^(-e)*c,16)),j=n+2)
			)
		);
		j++
	);
	if(j==n+1,return([1,t]))
);
for(i=0,k,
	X=2^i;
	for(Y=0,2^(k-i)-1,
		t=[X,Y];
		if(valuation([X,Y],2)==0,
			j=1;
			while(j<n+1 && !c=m[j,]*t~,j++);
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

