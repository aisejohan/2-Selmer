\\ Buhler-Gross algorithm for computing L-series derivatives

\\ Implementation by tom@womack.net
\\ Repaired and improved by jec@maths.nottingham.ac.uk

\\ This is practical for any rank and conductors <= 10^10 in 10
\\ minutes per L' calculation (on P3/500E).

default(primelimit,2000000);

\\ Section 0: prime-number utility which really should be
\\ built-in. This version due to Karim (reference for upper and lower
\\ bounds is Rosser & Schoenfeld, Approximate formulas for some functions
\\ of prime numbers, Illinois J Math 6 (1962) pp 64-94)

pith(x) =
{
  local(l,r,lx);
  if (x <= 2, return (x==2));
  lx = log(x);
  l = floor(x/lx);
  r = floor(l * (1 + 3/(2*lx)));
  while (r-l>1,
    m = (l+r)>>1;
    if (prime(m)<=x, l=m, r=m));
  l;
}

\\ Section 1: G_r

\\ We need a load of Magic Constants, which are the coefficients
\\ of the Taylor expansion of GAMMA(1+s). Now, [BGZ] we have
\\ log(Gamma(1+s)) = -gamma s + sum(n=2,infinity,(-1)^n*zeta(n)/n)
\\ and Pari can do polynomial exponential

default(seriesprecision,30);
ellld_ll = -Euler*s + sum(n=2,30,(-1)^n * zeta(n) / n * s^n);
ellld_mm = exp(ellld_ll);
ellld_pn_coeffs=vector(25,i,polcoeff(ellld_mm,i-1));

ellld_P(r,t)=local(m);sum(m=0,r,ellld_pn_coeffs[r-m+1]*t^m/factorial(m));

{ellld_G(r,x)=
   local(ss,oldss,n,sn);
   if (x>30,return(0));
   ss = ellld_P(r,log(1/x));
   oldss = ss-1;n=1;sn=-((-1)^r);
   while (abs(ss-oldss)>1e-15,
     oldss = ss;
     ss = ss + sn*x^n/(n^r*factorial(n));
     n=n+1;sn=-sn;
  );
  2*ss;
}

\\ Section 2: the BG recursion
\\ Some global variables are necessary, but I've named them
\\ to try to avoid contaminating the namespace *too* much.

ellld_getap(r)=if(r<=ellld_pnum,ellld_ap[r],ellap(est,prime(r)))
ellld_getp(r)=if(r<=ellld_pnum,ellld_p[r],prime(r))

{
BGadd(n,i,a,last_a) =
	local(j,j0,next_a);
	if(a==0,j0=i,ellld_sum=ellld_sum+a*ellld_G(ellld_r,ellld_X*n)/n;j0=1;);
	j=j0;
	while(j<=i && ellld_getp(j)*n <= ellld_bnd,
		next_a=a*ellld_getap(j);
		if(j==i && ellld_N%ellld_getp(j)!=0,
			next_a = next_a - ellld_getp(j)*last_a);
		BGadd(ellld_getp(j)*n,j,next_a,a);
		j=j+1;
	)
}

{
elllderiv(e,r,bnd) =
  local(cached_pmax,loop_p_num,ithresh,N,est);
  est	       = ellinit(e);
  cached_pmax  = sqrt(bnd);
  ellld_pnum   = pith(cached_pmax);
  loop_p_num   = pith(bnd);
  ellld_p      = primes(ellld_pnum);
  ellld_ap     = vector(ellld_pnum,i,ellap(est,ellld_p[i]));
\\ Version 2.0.19 has acknowledged bug in ellap for p=2
  ellld_ap[1]  = ellak(est,2);
  N	       = ellglobalred(e)[1];

  ellld_N      = N;
  ellld_X      = 2*Pi/sqrt(N);
  ellld_r      = r;
  ellld_bnd    = bnd;
  ellld_sum    = ellld_G(r,ellld_X);
  ithresh=2;
  for(i=1,loop_p_num,
\\    if(i==ithresh,
\\      print([i,ellld_sum]);
\\      ithresh=ithresh*2;
\\    );
    BGadd(ellld_getp(i),i,ellld_getap(i),1)
  );
  ellld_sum;
}

{
  ellanalyticrank(e)=
  local(rksgn,ZB,egrd,cond,localbnd,Lr1,RS,om);
  egrd=ellglobalred(e);
\\ We had to compute this anyway to get the conductor, and
\\ ellap gives wrong answers on non-minimal curves ...
  if(egrd[2]!=[1,0,0,0],
    print("INPUT CURVE NOT MINIMAL");
    e=ellchangecurve(e,egrd[2]);
    print("MINIMISED TO ",[e[1],e[2],e[3],e[4],e[5]]));
  cond=egrd[1];
\\ No, I can't justify this, but it should be big enough ...
  localbnd=round(2*sqrt(cond));
  print("Summing ",localbnd," a_n terms");
  rksgn=ellrootno(e);
  if(rksgn==-1,print("Rank is odd"),print("Rank is even"));  
  ZB=(1-rksgn)/2;
  Lr1=0;
\\ And I'm sure it's conjectural that this is close-enough to 0
  while(abs(Lr1)<0.001,
    Lr1=elllderiv(e,ZB,localbnd);
    print("L^(",ZB,")=",Lr1);
    ZB=ZB+2);
\\ Also compute the quantity reg*sha from the BSD conjecture
  om=e[15]; if(e.disc>0,om=2*om);
  [ZB-2,Lr1,(Lr1*elltors(e)[1]^2)/(om*egrd[3])];
}

{ellRegSha(e,r)=
  local(egrd,om,Lfn);
  om=e[15];if(e.disc>0,om=2*om);
  egrd=ellglobalred(e);
  Lfn=elllderiv(e,r,2*sqrt(egrd[1]));
  (Lfn*elltors(e)[1]^2)/(om*egrd[3]);
}  

ellSha(e,pts)=ellRegSha(e,length(pts))/matdet(ellheightmatrix(e,pts));
