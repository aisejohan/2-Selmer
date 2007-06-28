\r ell.gp

DEBUGLEVEL=0;

johan_lim=1000;
johan_a=1;
johan_b=2;
johan_got=0;
johan_sofar=vector(0);
johan_bnf=bnfinit(y);

johan_do_it_once()=
{
	if(johan_got>0,\
		johan_prime=nextprime(johan_sofar[johan_got]+1)\
	,\
		johan_prime=2\
	);
	johan_stop=0;
	johan_succes=0;
	johan_mask=0;
	while(!johan_stop,\
		johan_d=johan_prime;\
		for(johan_i=1,johan_got,\
			if(bittest(johan_mask,johan_i),\
				johan_d=johan_sofar[johan_i]*johan_d;\
			)\
		);\
		johan_rank=\
			complete(johan_bnf,johan_d*johan_a,johan_d*johan_b);\
		if(johan_rank,\
			if(johan_mask == 2^johan_got - 1,\
				print("For p=",johan_prime,\
					" the rank is: ",johan_rank)\
			);\
			print("Stop at:"\
				"\t johan_prime= ",johan_prime,\
				"\t johan_mask= ",johan_mask,\
				"\t johan_rank= ",johan_rank);\
			johan_prime=nextprime(johan_prime+1);\
			johan_mask=0;\
			if(johan_prime>johan_lim,johan_stop=1)\
		);\
		johan_mask++;\
		if(johan_mask == 2^johan_got,\
			johan_success=1;\
			johan_stop=1;\
		)\
	);
	if(johan_success,\
		johan_got++;\
		johan_sofar=concat(johan_sofar,[johan_prime])\
	);
	return(johan_success);
}
