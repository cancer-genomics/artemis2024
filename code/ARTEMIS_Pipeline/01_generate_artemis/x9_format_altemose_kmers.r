library("Biostrings")
library(tidyverse)
library(data.table)


fastaFile <- readDNAStringSet("x9_altemosereference24.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)

df$type<-sapply(str_split(df$seq_name,"_"),"[",1)

df$header<-paste0(">",df$seq_name)
f<-unique(df$type)
for (j in 1:length(f)) {
	d<-df %>% filter(type == f[j])
	d<-d %>% select(header,sequence)
	D <- do.call(rbind, lapply(seq(nrow(d)), function(i) t(d[i, ])))
	write.table(D, paste0("./final_kmers/",f[j],".fasta"),row.names = FALSE, col.names = FALSE, quote = FALSE)

}

