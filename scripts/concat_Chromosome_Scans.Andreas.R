library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
data_path = as.character(args[1])

# merge scan
files = dir(data_path, pattern = "R.pseudoscan.*.txt") # get file names
files = grep("chr",files,value=TRUE)
df1 = files %>%
  # read in all the files, appending the path before the filename
  map(~ as_tibble(read.table(file.path(data_path, .)))) %>% 
  reduce(rbind)

write.table(df1,paste0(data_path,"/R.pseudoscan.txt"))

# merge means
files = dir(data_path, pattern = "R.meansBySample.*.txt") # get file names
files = grep("chr",files,value=TRUE)

df2 = files %>%
  # read in all the files, appending the path before the filename
  map(~ as_tibble(read.table(file.path(data_path, .)))) %>% 
  reduce(rbind)

write.table(df2,paste0(data_path,"/R.meansBySample.txt"))


library(patchwork)

xx2 = df1 %>% mutate(chr2=recode(chr, "chrX" = "X", "chr2L" = "II", "chr2R" = "II", "chr3L" = "III", "chr3R" = "III"))
F1 = ggplot(xx2,aes(x=cM,y=Wald_log10p)) +
		geom_point(pch=15,cex=0.25) +
		facet_wrap(~ chr2, nrow=1)
		
F2 = ggplot(xx2 %>% mutate(Mb=pos/1e6),aes(x=Mb,y=Wald_log10p)) +
		geom_point(pch=15,cex=0.25) +
		facet_wrap(~ chr, nrow=1)
		
png(paste0(data_path,"/Ageing.png"), width=8, height=8, units="in", res=600)
F1 + F2 + plot_layout(ncol=1)
dev.off()

F3 = ggplot(xx2 %>% mutate(Mb=pos/1e6),aes(x=Mb,y=Wald_log10p)) +
		# geom_point(pch=15,cex=0.25) +
		# geom_smooth(se = FALSE, method = "gam", formula = y ~ s(x, bs = "cs", fx = TRUE, k = 200)) +
		geom_line() +
		facet_wrap(~ chr, ncol=1)

png(paste0(data_path,"/Ageing2.png"), width=8, height=10, units="in", res=600)
F3
dev.off()

F4 = ggplot(xx2 ,aes(x=cM,y=Wald_log10p)) +
                # geom_point(pch=15,cex=0.25) +
                # geom_smooth(se = FALSE, method = "gam", formula = y ~ s(x, bs = "cs", fx = TRUE, k = 200)) +
                geom_line() +
                facet_wrap(~ chr, scales = "free_x", ncol=1)

png(paste0(data_path,"/Ageing3.png"), width=8, height=10, units="in", res=600)
F4
dev.off()
