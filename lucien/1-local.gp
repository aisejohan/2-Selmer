/* Local computations. */

MAXPROB = 20;

/* Silly probabilistic attempt. */
/* The equations are:
	y1^2 = u*v*d*z
	y2^2 = u^(-1)*d*(z-a)
	y3^2 = v^(-1)*d*(z+b)
Can also write four quadratic equations
	y1^2 = u*v*d*z
	y2^2 = u^(-1)*d*(z-a*w)
	y3^2 = v^(-1)*d*(z+b*w)
	y4^2 = w
*/
has_points_at_p(a,b,d,u,v,p)=
{
	local(m,t,u);
	if(isprime(p),, error("Not a prime!"));
	if(valuation(u,p) > 1 || valuation(v,p) > 1 || valuation(d,p) > 1,
		error("reduce powers of p first."));
	m=matrix(4,2);
	m[1,1]=u*v*d;
	m[1,2]=0;
	m[2,1]=u^(-1)*d;
	m[2,2]=u^(-1)*d*(-a);
	m[3,1]=v^(-1)*d;
	m[3,2]=v^(-1)*d*b;
	m[4,1]=0;
	m[4,2]=1;
	return(recursive_has_points_at_p(m,p,5))
}

recursive_has_points_at_p(m,p,c)=
{
	local(t,i);
	c--;
	t=four_linear_forms(m,p);
	if(t,,
		i=1;
		while(i<5,
			t=[Mod(random(p),p),Mod(random(p),p)];
			i=1;
			while(i < 5 && issquare(m[i,1]*t[1]+m[i,2]*t[2]), i++)
		);
		t=lift(t)
	);
	if(Mod(t[1],p),
		u=[t[1],0;t[2],p]
	,
		u=[t[1],p;t[2],0]
	);
	m=m*u;
	return(1)
}

/* Answers the question: Can there be a solution, and if so:
what should it look like mod p? Coefficients in Q. */
four_linear_forms(m,p)=
{
	local(e,i,j,t);
	t=vector(2);
	for(i=1,4,
		e=valuation(m[i,],p);
		if(Mod(e,2),
			if(t && 
			(
			t[1]*Mod(-m[i,2]*p^(-e),p) -
			t[2]*Mod(m[i,1]*p^(-e),p)
			),
				return([0,[0,0]])
			,
				t=[Mod(-m[i,2]*p^(-e),p),Mod(m[i,1]*p^(-e),p)]
			)
		)
	);
	for(i=1,4,
		e=valuation(m[i,],p);
		m[i,]=p^(-e)*m[i,]
	);
	for(i=1,4,for(j=i+1,4,
		if(Mod(m[i,1]*m[j,2]-m[i,2]*m[j,1],p),,
			if(Mod(m[i,1],p),
				if(!issquare(Mod(m[j,1]/m[i,1],p)),
					return([-m[i,2]*p^(-e),m[i,1]*p^(-e)])
				)
			,
				if(!issquare(Mod(m[j,2]/m[i,2],p)),
					return([-m[i,2]*p^(-e),m[i,1]*p^(-e)])
				)

			)
		)
	));
	return([0,0])
}
	
