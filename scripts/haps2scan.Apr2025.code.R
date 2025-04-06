doscan2 = function(df,chr,Nfounders){
	sexlink = 1
	if(chr=="chrX"){ sexlink=0.75 }

	# I tested with xx2$data[[1]]
	df2 = df %>%
		unnest(cols = c(sample, Groups, Haps, Err, Names)) %>%
		left_join(design.df, join_by(sample==bam)) %>%
		filter(!is.na(TRT))
	
	# only analyze data for which all founders are discernable..
	allFounders = as.numeric(df2 %>% mutate(mm = max(unlist(Groups))) %>% summarize(max(mm)))	

	ll = list(Wald_log10p = NA, Pseu_log10p = NA, Falc_H2 = NA, Cutl_H2 = NA, avg.var = NA)
	if(allFounders!=Nfounders){ return(ll) }

	##  now cases where all founders are OK
	##  now collapse any pure replicates.  This is tidy ugly.  But I feel there 
	##  is value in keeping dataframe columns as lists...
	df3 = df2 %>%
		select(-Groups) %>%
		group_by(TRT,REP) %>%	
		summarise(Err_mean = list(reduce(map(Err, ~as.matrix(.x)), `+`)/length(Err)),
			Haps_mean = list(reduce(map(Haps, ~as.vector(.x)), `+`)/length(Haps)),
			Names = list(first(Names)),
			Num_mean = sexlink*mean(Num)) %>%
		rename(Haps=Haps_mean,Num=Num_mean,Err=Err_mean)

	## these summaries of the data are pretty useful for tests
	p1 = df3 %>% filter(TRT=="C") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t() 
	row.names(p1) <- NULL
	p2 = df3 %>% filter(TRT=="Z") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
	row.names(p2) <- NULL
	covar1 = do.call(abind, c(df3 %>% filter(TRT=="C") %>% pull(Err), along = 3))
	covar2 = do.call(abind, c(df3 %>% filter(TRT=="Z") %>% pull(Err), along = 3))
	nrepl = df3 %>% filter(TRT=="C") %>% nrow()
	nrepl == df3 %>% filter(TRT=="Z") %>% nrow()
	N1 = df3 %>% filter(TRT=="C") %>% pull(Num)
	N2 = df3 %>% filter(TRT=="Z") %>% pull(Num)

	wt=wald.test3(p1,p2,covar1,covar2,nrepl,N1,N2)
	Wald_log10p = -log10(wt$p.value)
	Pseu_log10p = pseudoN.test(p1,p2,covar1,covar2,nrepl,N1,N2)

	af_cutoff = 0.01     # 1% --- heritability estimators can be off for really low allele frequencies
	temp = Heritability(p1, p2, nrepl, ProportionSelect, af_cutoff)
	Falc_H2 = temp$Falconer_H2
	Cutl_H2 = temp$Cutler_H2

	ll = list(Wald_log10p = Wald_log10p, Pseu_log10p = Pseu_log10p,
			Falc_H2 = Falc_H2, Cutl_H2 = Cutl_H2, avg.var = wt$avg.var)
	ll
	}

xx1 = readRDS(filein)
Nfounders=length(xx1$Groups[[1]][[1]])
ProportionSelect = design.df %>% filter(TRT=="Z") %>% select(REP,Proportion) %>% arrange(REP)

bb1 = xx1 %>%
#	head(n=100) %>%
	group_by(CHROM,pos) %>%
	nest() %>%
	mutate(out = map2(data, CHROM, doscan2, Nfounders=Nfounders)) %>%
	unnest_wider(out)
bb2 = bb1 %>% select(-data) %>% rename(chr=CHROM)
bb3 = add_genetic(bb2)

# I drop the loci for which the scan gives and NA
bb4 = bb1 %>%
	filter(!is.na(Pseu_log10p)) %>%
	select(-c(Wald_log10p, Pseu_log10p, Falc_H2, Cutl_H2, avg.var, data)) %>%
	left_join(xx1) %>%
	select(-c(Err,Groups)) %>%
	unnest(c(sample,Haps,Names)) %>%
	unnest(c(Haps,Names)) %>%
	rename(chr=CHROM,pool=sample,freq=Haps,founder=Names) %>%
	left_join(design.df, by=c("pool"="bam")) %>%
	select(c(chr,pos,TRT,REP,REPrep,freq,founder)) %>%
	filter(!is.na(TRT)) %>%
	group_by(chr,pos,TRT,REP,founder) %>%
	summarize(freq=mean(freq,na.rm=TRUE))

write.table(bb3, fileout)
write.table(bb4, fileout_meansBySample)

