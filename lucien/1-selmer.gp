/* Some stuff about the (twisted) Frey curve:

	E_d : d*y^2 = x*(x-a)*(x+b)

Let E^d be the curve given by

	E^d : y^2 = x*(x-d*a)*(x+d*b)

There is a morphism

	E_d ---> E^d
	(x,y) \mapsto (d*x,d^2*y).

Check: */
tmp = subst(subst(y^2-x*(x-d*a)*(x-d*b),x,d*x),y,d^2*y) - \
d^3*(d*y^2-x*(x-a)*(x-b));
if(tmp,error("Quel horreur!"))

/* So we can use the Weierstrass form E^d to compute invariants
of the twisted curve E_d. Write the equation of E^d as 

	y^2 = f(x).

*/

f=x*(x-d*a)*(x+d*b);

/* We use eee=E^d in pari/gp. */
eee=ellinit([0,polcoeff(f,2,x),0,polcoeff(f,1,x),polcoeff(f,0,x)],1);

/* The two torsion points on E^d. */
e1=[0,0];
e2=[d*a,0];
e3=[-d*b,0];

/* Generic point on E^d. */
generic=[x,y];

/* The x-coordinates of adding the two torsion points on E_d. */
/* First add on E^d.*/
tmp=elladd(eee,generic,e1)[1];
x1=subst(tmp,y,0)+polcoeff(tmp,2,y)*f;
/* Rework in terms of x-coordinate on E_d. */
x1=subst(x1,x,d*x)/d;

tmp=elladd(eee,generic,e2)[1];
x2=subst(tmp,y,0)+polcoeff(tmp,2,y)*f;
x2=subst(x2,x,d*x)/d;

tmp=elladd(eee,generic,e3)[1];
x3=subst(tmp,y,0)+polcoeff(tmp,2,y)*f;
x3=subst(x3,x,d*x)/d;

/* The y-coordinates of adding the two torsion points on E_d. */
/* First add on E^d.*/
tmp=elladd(eee,generic,e1)[2];
y1=y*(polcoeff(tmp,1,y)+polcoeff(tmp,3,y)*f);
/* Rework in terms of x-coordinate on E_d. */
y1=subst(subst(y1,x,d*x),y,d^2*y)/d^2;

tmp=elladd(eee,generic,e2)[2];
y2=y*(polcoeff(tmp,1,y)+polcoeff(tmp,3,y)*f);
y2=subst(subst(y2,x,d*x),y,d^2*y)/d^2;

tmp=elladd(eee,generic,e3)[2];
y3=y*(polcoeff(tmp,1,y)+polcoeff(tmp,3,y)*f);
y3=subst(subst(y3,x,d*x),y,d^2*y)/d^2;

/* Check that the maps satisfy the composition law. */
if(x1 - subst(x2,x,x3),error("Quel horreur!"))

/* The x-coordinate of multiplication by 2 on E_d. */
/* Again first compute on E^d and then rework. */
tmp=elladd(eee,generic,generic)[1];
varphi=polcoeff(tmp,0,y)+polcoeff(tmp,-2,y)*f^(-1);
varphi=subst(varphi,x,d*x)/d;

/* The y-coordinate of multiplication by 2 on E_d. */
/* Again first compute on E^d and then rework. */
tmp=elladd(eee,generic,generic)[2];
varpsi=polcoeff(tmp,1,y)*f^2*y^(-3)+\
	polcoeff(tmp,-1,y)*f*y^(-3)+\
	polcoeff(tmp,-3,y)*y^(-3);
varpsi=subst(subst(varpsi,x,d*x),y,d^2*y)/(d^2);

/* Check that varphi is compatible with adding 2-torsion points. */
if(varphi - subst(varphi,x,x1),error("Quel horreur!"))
if(varphi - subst(varphi,x,x2),error("Quel horreur!"))
if(varphi - subst(varphi,x,x3),error("Quel horreur!"))

/* Now we can express varphi in terms of x+xi. */
g1 = (u^2 + 4*a*b)/(4*(u-a+b));
if(varphi - subst(g1,u,x+x1), error("Quel horreur!"))

g2 = u^2/(4*(u+b));
if(varphi - subst(g2,u,x+x2), error("Quel horreur!"))

g3 = u^2/(4*(u-a));
if(varphi - subst(g3,u,x+x3), error("Quel horreur!"))

/* At this point, for any z in P^1, the torsor for E_d[2] corrsponding
to z is simply given by the equation

	varphi = z.

This length 4 scheme has an action of the maps xi i=1,2,3. The torsor
is trivial iff the equations

	gi(u) = z

have solutions in u for i=1,2,3. This is a quadratic equation in each case.
We compute their disciminants. */

d1 = 4*(z-a)*(z+b);
if(4*(u-a+b)*(g1 - z) - ((u-2*z)^2 - d1), print("Quel horreur!"))

d2 = 4*z*(z+b);
if(4*(u+b)*(g2 - z) - ((u-2*z)^2 - d2), print("Quel horreur!"))

d3 = 4*z*(z-a);
if(4*(u-a)*(g3 - z) - ((u-2*z)^2 - d3), print("Quel horreur!"))

/* If z is the x-coordinate of a point on the curve E_d the
condition that d1 be a square is equivalent to either of the
conditions

	(z-a)*(z+b) is a square, or
	d*z is a square.

This is where the relevance of d comes in. Similar for d2, d3. 
Consider the curve C_d defined by the equations

	y1^2 = d*z
	y2^2 = d*(z-a)
	y3^2 = d*(z+b)

This is a genus one curve and it maps to E_d by the map

	C_d ---> E_d
	(z, y1, y2, y3) \mapsto (z, y1*y2*y3/d^2).

Of course, by the above, there is an isomorphism

	ISO : E_d ---> C_d
	(x, y) \mapsto (varphi, 
			(x^2+a*b)/(2*y), 
			(x^2-2*a*x-a*b)/(2*y),
			(x^2+2*b*x-a*b)/(2*y))

such that the composition of the maps is the multiplication by 2 map.
Check: */
ISO=[varphi,(x^2+a*b)/(2*y), (x^2-2*a*x-a*b)/(2*y), (x^2+2*b*x-a*b)/(2*y)];
tmp = d*ISO[1] - ISO[2]^2;
tmp = polcoeff(tmp,0,y) + polcoeff(tmp,-2,y)*(x*(x-a)*(x+b)/d)^(-1);
if(tmp,print("Quel horreur!"))
tmp = d*(ISO[1]-a) - ISO[3]^2;
tmp = polcoeff(tmp,0,y) + polcoeff(tmp,-2,y)*(x*(x-a)*(x+b)/d)^(-1);
if(tmp,print("Quel horreur!"))
tmp = d*(ISO[1]+b) - ISO[4]^2;
tmp = polcoeff(tmp,0,y) + polcoeff(tmp,-2,y)*(x*(x-a)*(x+b)/d)^(-1);
if(tmp,print("Quel horreur!"))
if(varpsi - ISO[2]*ISO[3]*ISO[4]/d^2,print("Quel horreur!"))

/*
By the computations above, the point (x,y) is in the image of the
map C_d --> E_d if and only if it is divisble by 2 in E(K). In other
words, C_d --> E_d is a model for the multiplication by 2 map. The
action of {1,-1} x E[2] is given by the formulae:

	[-1]	turns into	yi \mapsto -yi (i=1,2,3)
	x1	turns into	yi \mapsto -yi (i=2,3)
	x2	turns into	yi \mapsto -yi (i=1,3)
	x3	turns into	yi \mapsto -yi (i=1,2)

Check: */
if(subst(subst(ISO[2],x,x1),y,y1) - ISO[2],error(""))
if(subst(subst(ISO[3],x,x1),y,y1) + ISO[3],error(""))
if(subst(subst(ISO[4],x,x1),y,y1) + ISO[4],error(""))
if(subst(subst(ISO[2],x,x2),y,y2) + ISO[2],error(""))
if(subst(subst(ISO[3],x,x2),y,y2) - ISO[3],error(""))
if(subst(subst(ISO[4],x,x2),y,y2) + ISO[4],error(""))
if(subst(subst(ISO[2],x,x3),y,y3) + ISO[2],error(""))
if(subst(subst(ISO[3],x,x3),y,y3) + ISO[3],error(""))
if(subst(subst(ISO[4],x,x3),y,y3) - ISO[4],error(""))

/* For a fixed d we can twist by cocyles in 

	H^1(K, E_d[2]) \cong  (K^* x K^*)/2

by letting (u,v) correspond to the twist C_{d,u,v}

	y1^2 = u*v*d*z
	y2^2 = u^(-1)*d*(z-a)
	y3^2 = v^(-1)*d*(z+b)

of the curve C_d with twisted morphism

	C_{d,u,v} ---> E_d
	(z, y1, y2, y3) \mapsto (z, y1*y2*y3/d^2).

*/

/* Trivial elements:

	(0, 0, 1, 1)	 	on C_{d,-a*d,b*d} maps to (0,0)
	(a, d*a*(a+b), 0, 1)	on C_{d,a*(a+b),(a+b)*d} maps to (a,0)
	(-b, d*b*(a+b), 1, 0)	on C_{d,-(a+b)*d,b*(a+b)} maps to (-b,0)

*/

/* The curve C_{d,u,v} maps to the rational curve R_{u,v}

	(u*v)^(-1)*(a+b)*y1^2 = u*b*y2^2 + v*a*y3^2.

This is equivalent to the equation

	y1^2 = u^2*v*(a+b)*b*y2^2 + u*v^2*a*(a+b)*y3^2.

The quadratic form associated to the Quaternion algebra 
K<i,j>/(i^2-f,j^2-g,ij+ji) is fX^2+gY^2-fgZ^2 which is 
equivalent to (1/f)X^2+(1/g)Y^2-Z^2 and fX^2+gY^2-Z^2. Hence the
symbol associated the the rational curve R_{u,v} is

	( v(a+b)b, ua(a+b) )_2.

The function field of this rational curve is the field

	K(y2/y1)[y3/y1]/( u*b*(y2/y1)^2 + v*a*(y3/y1)^2 - (a+b)/(u*v) ).

*/

/*
Solving the equations

	y1^2 = u*v*d*z
	y2^2 = u^(-1)*d*(z-a)
	y3^2 = v^(-1)*d*(z+b)

is equivalent to solving the system of homogeneous equations

	Y1^2 = u*v*d*Z
	Y2^2 = u^(-1)*d*(Z-a*W)
	Y3^2 = v^(-1)*d*(Z+b*W)
	Y4^2 = W

and this in turn is equivalent to solving the pair of 
quadratic equations

	u*Y2^2 - v*Y3^2 = (-a*d - b*d)*W = -(a+b)*d*Y4^2
	b*u*Y2^2 + a*v*Y3^2 = (b*d + a*d)*Z = (a+b)*u^(-1)*v^(-1)*Y1^2

or

	                          u*Y2^2 - v*Y3^2 + (a+b)*d*Y4^2 = 0
	-(a+b)*(uv)^(-1)*Y1^2 + b*u*Y2^2 + a*v*Y3^2  = 0

Given a solution (Y1,Y2,Y3,Y4) we get:

	W  = Y4^2
	Z  = Y1^2/(u*v*d)
	y1 = Y1/Y4
	y2 = Y2/Y4
	y3 = Y3/Y4
	z  = Z/W
	x  = z
	y  = y1*y2*y3/d^2

*/
