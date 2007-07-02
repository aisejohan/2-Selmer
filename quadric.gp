X;
Y;
Z;

equations(G,P)=
{
local(i,j,v,A);
P=P/content(P);
i=1;
A=abs(P[1]);
j=2;
while(j<5,
	if(abs(P[j]) > A,
		A=abs(P[j]);
		i=j
	);
	j++
);
v=vector(4);
if(i == 1, v=[0,X,Y,Z]);
if(i == 2, v=[X,0,Y,Z]);
if(i == 3, v=[X,Y,0,Z]);
if(i == 4, v=[X,Y,Z,0]);
return((v*G*v~)*P-2*(P*G*v~)*v)
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

find_smaller(G,P,lijst,H=0)=
{
local(t,X,first);
P=P/content(P);
if(H == 0, 
	H=sum(i=1,4,P[i]^2);
	first=1
,
	first=0
);
for(i=1,matsize(lijst)[2],
	X=lijst[i];
	t=-2*X*G*P~/(X*G*X~);
	X=P+t*X;
	X=X/content(X);
	if(H > sum(j=1,4,X[j]^2),
		return(X)
	);
	if(first,
		find_smaller(G,P,lijst,H)
	)
);
return(0)
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
	M[,3]=lattice[,3];
	if(matrank(M) == 3,
		v1=lattice[,1];
		v2=lattice[,3]
	,
		M[,2]=lattice[,2];
		if(matrank(M) == 3,
			v1=lattice[,2];
			v2=lattice[,3]
		,
			error("Should not happen!")
		)
	)
);
Q=matrix(2,2);
Q[1,1]=(v1~)*G*v1;
Q[1,2]=(v1~)*G*v2;
Q[2,1]=(v2~)*G*v1;
Q[2,2]=(v2~)*G*v2;
return(Q)
}

