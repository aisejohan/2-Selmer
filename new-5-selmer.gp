\r weak_selmer.gp.run

global(j_a,j_b,j_lim,j_got,j_sofar,j_sofar_len,j_sofar_bounds);

do_it_all()=
{
j_a=1;
j_b=2;
j_lim=5000;
j_got=0;
j_sofar=vector(0);
j_sofar_len=vector(0);
j_sofar_bounds=vector(0);
j_setup(3);
j_do_it_once(5);
j_do_it_once(7);
j_do_it_once(9);
return(0)
}

j_compute(j_d)=
{
	local(j_return);
	j_return=global_linear_conditions(j_a,j_b,j_d);
	if(j_return > 10,
		print("Here the 2-Selmer rank is: ",j_return);
		write("BIGSELMER",
		"a=",j_a,"; b=",j_b,"; d="j_d,"; selmer=",j_return,";")
	);
	return(j_return)
}

j_setup(j_bound)=
{
	local(j_prime,j_rank);
	j_prime=2;
	j_sofar=concat(j_sofar,[vector(0)]);
	j_sofar_len=concat(j_sofar_len,[0]);
	j_sofar_bounds=concat(j_sofar_bounds,[j_bound]);
	while(j_prime<j_lim,
		j_rank=j_compute(j_prime);
		if(j_rank > j_bound,
			print(j_rank);
			j_sofar_len[1]++;
			j_sofar[1]=concat(j_sofar[1],[j_prime])
		);
		j_prime=nextprime(j_prime+1)
	);
	j_got=1
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

j_do_it_once(j_bound)=
{
	local(j_tmp,j_factors,j_j,j_prime,j_ok,j_n,j_spot,j_rank);
	j_sofar=concat(j_sofar,[vector(0)]);
	j_sofar_len=concat(j_sofar_len,[0]);
	j_sofar_bounds=concat(j_sofar_bounds,[j_bound]);
	for(j_i=1,j_sofar_len[j_got],
print("Doing ",j_i," out of ",j_sofar_len[j_got]);
		j_tmp=j_sofar[j_got][j_i];
		j_factors=factor(j_tmp)[,1];
		j_j=find_spot(
			j_sofar[1],
			j_factors[j_got]+1);
		while(j_j <= j_sofar_len[1],
			j_prime=j_sofar[1][j_j];
			j_ok=1;
			for(j_k=1,j_got,
				j_n=j_tmp*j_prime/j_factors[j_k];
				j_spot=find_spot(j_sofar[j_got],j_n);
				if(j_spot > j_sofar_len[j_got] || 
					j_n != j_sofar[j_got][j_spot],
					j_ok=0;
					j_k=j_got;
				)
			);
			if(j_ok,
				j_rank=j_compute(j_tmp*j_prime);
				if(j_rank > j_bound,
					j_sofar_len[j_got+1]++;
					j_sofar[j_got+1]=
						concat(j_sofar[j_got+1],
							[j_tmp*j_prime])
				)
			);
			j_j++
		)	
	);
	j_sofar[j_got+1]=vecsort(j_sofar[j_got+1]);
	j_got++
}
