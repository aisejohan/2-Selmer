\r weak_selmer.gp.run

j_compute(j_a,j_b,j_d)=
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

j_setup(j_a,j_b,j_bound,j_lim)=
{
	local(l,len,j_prime,j_rank);
	l=[];
	len=0;
	j_rank=j_compute(j_a,j_b,-1);
	if(j_rank > j_bound,
		print("For -1 we get: ",j_rank);
		len++;
		l=concat(l,[-1])
	);
	j_prime=2;
	while(j_prime<j_lim,
		j_rank=j_compute(j_a,j_b,j_prime);
		if(j_rank > j_bound,
			print(j_rank);
			len++;
			l=concat(l,[j_prime])
		);
		j_prime=nextprime(j_prime+1)
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
				if(j_rank > j_bound,
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
