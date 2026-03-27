library(seqinr)
library(ggplot2)

#DIR <- '/mnt/trcanmed/snaketree/prj/Sanger_traces/local/share/data/run_0326'
fsa_file <- snakemake@input[['fsa']]
outplot <- snakemake@output[['plot']]
n_chan <- as.numeric(snakemake@params[['n_channels']])
log_f <- snakemake@log[['log']]
### Code for manual inspection of fsa files to determine which is the ladder channel for each of them 
#d_fsas <- read_fsa(fsas) # this is a trace function that under the hood uses read.abif from the seqinr package,
# also read_fsa from Fragman uses it 
sink(log_f)
d0 <- read.abif(fsa_file, verbose=TRUE) # same as d_fsas[[1]]$fsa
sink()

intensities <- c()
channels <- c()
n <- 0
begin <- grep('DATA.1', names(d0$Data), fixed=TRUE)
for (i in seq(0, n_chan-1)) {
  n <- length(d0$Data[[i+begin]])
  intensities <- c(intensities, d0$Data[[i+begin]])
  channels <- c(channels, c(rep(paste0('channel_', i+1), n)))
}
pd <- data.frame(intensity=intensities, channel=channels, x=rep(seq(1, n), n_chan))

p1 <- ggplot(data=pd, aes(y=intensity, x=x, color=channel))+geom_line()+theme_bw(base_size=10)
ggsave(p1, file=outplot)


# Note that there are more d0$Data$DATA.1-4 (also DATA.5-8 in the first example) but only the first four with 'sensible' numbers

# These are directly the raw data stored in the files. Which fluorofore is there / where is the ladder is not clear
# it was in the fourth channel in the first data, I would check looking at the data directly where is the ladder.
# We need a ladder for each fsa file to map these signals (all channels of a single fsa file share the length
# of the signal in different channels) to basepairs
# 
