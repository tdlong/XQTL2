xx1 = readRDS(filein)
recodeTable = as_tibble(data.frame(sample=names_in_bam,pool=samples))
Nfounders=length(xx1$Groups[[1]][[1]])

bb1 = xx1 %>%
#	head(n=100) %>%
	group_by(CHROM,pos) %>%
	nest() %>%
	mutate(out = map2(data, CHROM, doscan, Nfounders=Nfounders)) %>%
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
	rename(chr=CHROM,pool=sample,freq=Haps,founder=Names)

write.table(bb3,fileout)
write.table(bb4, fileout_meansBySample)

