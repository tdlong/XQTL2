window_size_multiplier = function(chr,position){
  # Convert position from bp to Mb
  pos_Mb <- position / 1e6
  
  # Maximum RR value (you'll need to fill this in for each chromosome)
  if(chr == "chr2L") {
    max_RR <- 4.135268  # the value you just calculated
    
    # Your 8th degree polynomial coefficients
    coeffs <- c( 0.7761563796515, 0.93380007449519, 0.433427661585053, -0.224762137968782, 0.0387093059807801, -0.00335077841277375, 0.000156654074237256, -3.78015959460574e-06, 3.70105691705139e-08 )
    
    # Calculate RR at the given position
    RR <- coeffs[1] + coeffs[2] * pos_Mb + coeffs[3] * pos_Mb^2 + 
          coeffs[4] * pos_Mb^3 + coeffs[5] * pos_Mb^4 + coeffs[6] * pos_Mb^5 + 
          coeffs[7] * pos_Mb^6 + coeffs[8] * pos_Mb^7 + coeffs[9] * pos_Mb^8
    
    return(max_RR / max(0.1, RR))
    
  } else if(chr == "chr2R") {
    max_RR <- 3.760792  # the value you just calculated
    
    # Your 8th degree polynomial coefficients
    coeffs <- c( -0.628260149363033, 0.936354223696003, -0.393500382518582, 0.0569385095543652, -0.000780601646806972, -0.000396750443499961, 3.24836150059996e-05, -1.0000001368976e-06, 1.105750119529e-08 )

    # Calculate RR at the given position
    RR <- coeffs[1] + coeffs[2] * pos_Mb + coeffs[3] * pos_Mb^2 + 
          coeffs[4] * pos_Mb^3 + coeffs[5] * pos_Mb^4 + coeffs[6] * pos_Mb^5 + 
          coeffs[7] * pos_Mb^6 + coeffs[8] * pos_Mb^7 + coeffs[9] * pos_Mb^8
    
    return(max_RR / max(0.1, RR))
    
  } else if(chr == "chr3L") {
    max_RR <- 2.950546  # the value you just calculated
    
    # Your 8th degree polynomial coefficients
    coeffs <- c( 0.48384469066167, 2.74307261058466, -1.30341119964166, 0.293155520677805, -0.035039985867133, 0.00236275596110228, -9.05660725451677e-05, 1.84321638605809e-06, -1.54703865205797e-08 )

    # Calculate RR at the given position
    RR <- coeffs[1] + coeffs[2] * pos_Mb + coeffs[3] * pos_Mb^2 + 
          coeffs[4] * pos_Mb^3 + coeffs[5] * pos_Mb^4 + coeffs[6] * pos_Mb^5 + 
          coeffs[7] * pos_Mb^6 + coeffs[8] * pos_Mb^7 + coeffs[9] * pos_Mb^8
    
    return(max_RR / max(0.1, RR))
    
  } else if(chr == "chr3R") {
     max_RR <- 3.413854  # the value you just calculated
    
    # Your 8th degree polynomial coefficients
    coeffs <- c( 0.0801229339474173, 0.0710663983833287, -0.134221170764257, 0.0474648132529972, -0.00692300110208564, 0.000535921800922228, -2.27159227968791e-05, 4.94302625176729e-07, -4.30809909384175e-09 )

    # Calculate RR at the given position
    RR <- coeffs[1] + coeffs[2] * pos_Mb + coeffs[3] * pos_Mb^2 + 
          coeffs[4] * pos_Mb^3 + coeffs[5] * pos_Mb^4 + coeffs[6] * pos_Mb^5 + 
          coeffs[7] * pos_Mb^6 + coeffs[8] * pos_Mb^7 + coeffs[9] * pos_Mb^8
    
    return(max_RR / max(0.1, RR))
  } else if(chr == "chrX"){
    max_RR <- 4.777405  # the value you just calculated
    
    # Your 8th degree polynomial coefficients
    coeffs <- c( 1.33605614233481, -4.72044955430154, 3.93867129150768, -1.08500846892869, 0.148155421160987, -0.0112461066139873, 0.000483366555582591, -1.10204966095983e-05, 1.03635540161062e-07 )
    
    # Calculate RR at the given position
    RR <- coeffs[1] + coeffs[2] * pos_Mb + coeffs[3] * pos_Mb^2 + 
          coeffs[4] * pos_Mb^3 + coeffs[5] * pos_Mb^4 + coeffs[6] * pos_Mb^5 + 
          coeffs[7] * pos_Mb^6 + coeffs[8] * pos_Mb^7 + coeffs[9] * pos_Mb^8
    
    return(max_RR / max(0.1, RR))
  } else {
    stop("Chromosome not recognized")
  }
}


nrow_subset = function(spotsdf, df3){
        df3 %>%
                filter(CHROM==spotsdf$CHROM &
                        POS > spotsdf$start &
                        POS < spotsdf$end &
                        (name %in% founders | name %in% names_in_bam)) %>%
                nrow()
        }




est_hap = function(spotsdf, df3){
        # spotsdf = spots$data[[1]]  testing
        temp_mat = df3 %>%
                filter(CHROM==spotsdf$CHROM &
                        POS > spotsdf$start &
                        POS < spotsdf$end &
                        (name %in% founders | name %in% names_in_bam)) %>%
                select(-c(CHROM,N)) %>%
                pivot_wider(names_from=name, values_from=freq) %>%
                pivot_longer(!c("POS",matches(founders)),names_to = "sample", values_to = "freq") %>%
                select(-POS)

        sample_mat = temp_mat %>%
                group_by(sample) %>%
                nest() %>%
                mutate(haps=map(data,est_hap2)) %>%
                select(-data) %>%
                unnest_wider(haps)

        sample_mat
        }
        
# function for estimating the haplotype for a given sample, called from est_hap
# returns Groups from cuttree, hap freq estimates, errors by founder
est_hap2 = function(sampdf){

        # sampdf = sample_mat$data[[1]]  # testing
        founder_mat = sampdf %>% select(matches(founders))
        Y = sampdf$freq
        good = !is.na(Y)
        Y = Y[good]
        founder_mat = founder_mat[good,] 
        m_founder_mat = as.matrix(founder_mat)
        Groups = cutree(hclust(dist(t(m_founder_mat))),h=h_cutoff)
        d = ncol(m_founder_mat)         
        out = lsei(A=m_founder_mat,B=Y, E=t(matrix(rep(1,d))),F=1,G=diag(d),H=matrix(rep(0.0003,d)),verbose=TRUE,fulloutput=TRUE)
        Haps = out$X
        Err = out$cov
        list(Groups=Groups,Haps=Haps,Err=Err,Names=names(Haps))
        }

# function for averaging error matrices across a sliding window
# takes a list of error matrix lists and returns averaged matrices
average_error_matrices = function(err_list){
        # err_list is a list where each element is Err[[i]] (a list of matrices, one per sample)
        # Get the number of samples from the first element
        n_samples = length(err_list[[1]])
        
        # Initialize result list
        averaged = vector("list", n_samples)
        
        # For each sample, average the corresponding matrices
        for(j in 1:n_samples){
                # Extract matrices for sample j across all rows in the window
                matrices = lapply(err_list, function(x) x[[j]])
                # Average the matrices
                averaged[[j]] = Reduce("+", matrices) / length(matrices)
        }
        
        return(averaged)
        }

df = lazy_dt(read.table(filein,header=TRUE))
df2 = df %>%
	pivot_longer(c(-CHROM,-POS), names_to = "lab", values_to = "count") %>%
	mutate(RefAlt = str_sub(lab,1,3)) %>%
	mutate(name = str_sub(lab,5)) %>%
	select(-lab) %>%
#	separate(lab, c("RefAlt", "name"), "_", extra = "merge") %>%
	pivot_wider(names_from = RefAlt, values_from = count) %>%
	mutate(freq = REF/(REF+ALT), N = REF+ALT) %>%
	select(-c("REF","ALT")) %>%
	as_tibble()

rm(df)
cat("df2 is now made\n")

# identify SNPs that are NOT problematic in the set of founders
good_SNPs = df2 %>%
	filter(name %in% founders) %>%
	group_by(CHROM,POS) %>%
	summarize(zeros=sum(N==0),notfixed=sum(N!=0 & freq > 0.03 & freq < 0.97),informative=(sum(freq)>0.05 | sum(freq) < 0.95)) %>%
	ungroup() %>%
	filter(zeros==0 & notfixed==0 & informative=="TRUE") %>%
	select(c(CHROM,POS))

# now subset the entire dataset for the good SNPs only
df3 = good_SNPs %>%
	left_join(df2, multiple = "all")

rm(df2)
cat("df3 is now made\n")
saveRDS(df3, file = rdsfile)

# df3 = readRDS(rdsfile)
# spots are the locations at which we will estimate haplotypes
# every <step> bp (i.e., 10kb) on the step (i.e, 0, 10, 20, ... kb)
# I define a window +/- size on those steps, and fix the ends
minpos = min(df3$POS)
maxpos = max(df3$POS)
myseq = seq(0,maxpos,step)
# Calculate window size multiplier for each position and create windows with variable sizes
spots = data.frame(CHROM=rep(mychr,length(myseq)), pos=myseq) %>%
  rowwise() %>%
  mutate(
    multiplier = window_size_multiplier(CHROM, pos),
    window_size = size * multiplier,
    start = pos - window_size,
    end = pos + window_size
  ) %>%
  ungroup() %>%
  # Filter windows: keep if pos is within chromosome, and if extending past ends,
  # keep only if shorter side is at least half the longer side
  mutate(
    keep_window = pos >= minpos & pos <= maxpos &  # midpoint must be within chromosome
      if_else(start < minpos, 
              (pos - minpos) >= 0.5 * (end - pos),  # left side: shorter >= 0.5 * longer
              TRUE) &
      if_else(end > maxpos,
              (maxpos - pos) >= 0.5 * (pos - start),  # right side: shorter >= 0.5 * longer
              TRUE)
  ) %>%
  filter(keep_window) %>%
  select(CHROM, pos, start, end)
# get rid of windows with fewer than 50 SNPs
# i.e., <50 SNPs in 100kb is pretty strange
UU = unique(df3$POS)
spots = spots %>%
                rowwise() %>%
                mutate(NN = sum(start < UU) - sum(end < UU)) %>%
                filter(NN >= 50) %>%
                select(-NN)

# this is the actual scan
# I guess this could also be slow...
# as it runs for L loci X S samples
spots2 = spots %>%
        group_nest(row_number()) %>%
        mutate(out = map(data,est_hap,df3)) %>%
		unnest(data) %>%
        select(-c(start,end)) %>%
        select(-`row_number()`) %>%
        unnest_wider(out) %>%
        # Add SWErr: sliding window averaged error matrices (±5 rows, 11 total)
        {
                n_rows = nrow(.)
                err_col = .$Err
                mutate(., 
                        SWErr = map(1:n_rows, function(i){
                                # Get window indices with boundary handling
                                start_idx = max(1, i - 10)
                                end_idx = min(n_rows, i + 10)
                                # Extract Err values for the window
                                err_window = err_col[start_idx:end_idx]
                                # Average the error matrices
                                average_error_matrices(err_window)
                        })
                )
        } %>%
        # Rename columns: Err -> nonSWErr, SWErr -> Err
        rename(nonSWErr = Err, Err = SWErr)

saveRDS(spots2,file = fileout)


