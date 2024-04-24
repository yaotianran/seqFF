#!/usr/bin/env Rscript
# seqFF.R [SE|PE] [bam_file]
options(tidyverse.quiet = TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressMessages(library(stringr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(tools))
suppressMessages(library(magrittr))
options(warn = -1)


args = commandArgs(trailingOnly = TRUE)
MODE = args[1]
BAM = args[2]


SCRIPT_FOLDER = dirname(this_file)



# = = = = = = = = = = = = = = = = = count = = = = = = = = = = = = = = = = =
TEMP_FOLDER = './'

this_file <- commandArgs() %>% tibble::enframe(name = NULL) %>% tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>% dplyr::filter(key == "--file") %>% dplyr::pull(value)
SCRIPT_FOLDER = paste0(dirname(this_file), '/')
CHROMS = paste0('chr', c(1:22, 'X', 'Y'))
r = load(paste0(DATA_FOLDER, 'seqFF.RData'))
bininfo.df = read_tsv(paste0(DATA_FOLDER, 'seqFF.tsv'), show_col_types = F, progress = F)  # 61927
temp_name = paste0('tmp.', as.character(sample(10000000000:99999999999, 1)))
file_name = str_split(basename(BAM), '\\.', simplify = T)[1, 1]

count_bam <- function(bam_file) {
   temp_folder = paste0(TEMP_FOLDER, '/', temp_name)

   if (!dir.exists(temp_name)) {
      dir.create(temp_name)
   }

   if (MODE == 'PE') {
      command_str = sprintf('samtools view -f 3 -F 3852 -o %s/%s.sam %s', temp_name, temp_name, bam_file)
   } else {
      command_str = sprintf('samtools view -F 2308 -o %s/%s.sam %s', temp_name, temp_name, bam_file)
   }

   message('counting ', BAM)
   r = system(command_str, intern = F)

   sam.file = sprintf('%s/%s.sam', temp_name, temp_name)
   sam.df = read_tsv(sam.file, col_names = F, quote = '\\', show_col_types = F, progress = F, col_types = cols())

   chrom.str = names(sort(table(sam.df$X3), decreasing = T))[1]
   if (!str_starts(chrom.str, 'chr')) {
      sam.df %<>% mutate(X3 = paste0('chr', X3))
   }
   sam.df %<>% dplyr::filter(X3 %in% CHROMS)

   sam.gr = makeGRangesFromDataFrame(sam.df, seqnames.field = 'X3', start.field = 'X4', end.field = 'X4')
   bininfo.gr = makeGRangesFromDataFrame(bininfo.df, seqnames.field = 'CHR', start.field = 'start', end.field = 'end')

   count.c = countOverlaps(bininfo.gr, sam.gr)
   bininfo.df %<>% mutate(counts = count.c)

   if (dir.exists(temp_folder)) {
      unlink(temp_folder, recursive = T)
   }
   rm(sam.df)
   gc()

   return(bininfo.df)
}

# WRSC predict
ff.pred <- function(gc.norm.bc.61927, B, mu, parameter.1, parameter.2){
   gc.norm.bc.61927[is.na(gc.norm.bc.61927)] <- 0
   gc.norm.bc <- gc.norm.bc.61927[grepl('chr[0-9]', names(gc.norm.bc.61927))]
   gc.norm.bc.c <- gc.norm.bc - mu
   y.hat <- matrix(c(1, gc.norm.bc.c), nrow = 1) %*% B
   y.hat.rep <- sum(y.hat, na.rm = T) / sum(gc.norm.bc)
   ff.hat <- (y.hat.rep + parameter.1) * parameter.2
   return(ff.hat)
}


# = = = = = = = = = = = = = = = = = process = = = = = = = = = = = = = = = = =
bininfo.df = count_bam(BAM)

autosomebins.bool <- bininfo.df$BinFilterFlag == 1 & bininfo.df$CHR != "chrX" & bininfo.df$CHR != "chrY"  # 50034
alluseablebins.bool <- bininfo.df$BinFilterFlag == 1  # 52611
autoscaledtemp <- bininfo.df$counts[autosomebins.bool] / sum(bininfo.df$counts[autosomebins.bool], na.rm = T)
allscaledtemp <- bininfo.df$counts[alluseablebins.bool] / sum(bininfo.df$counts[alluseablebins.bool], na.rm = T)

# additive loess correction
mediancountpergc.c <- tapply(autoscaledtemp, bininfo.df$GC[autosomebins.bool], function(x) median(x, na.rm = T))
loess.fitted.c <- predict(loess(mediancountpergc.c ~ as.numeric(names(mediancountpergc.c))), bininfo.df$GC[alluseablebins.bool])
normalizedbincount.c <- allscaledtemp + (median(autoscaledtemp, na.rm = T) - loess.fitted.c) # 52611
bincounts.c = rep(0, nrow(bininfo.df))
names(bincounts.c) = bininfo.df$binName
bincounts.c[alluseablebins.bool] <- (normalizedbincount.c / sum(normalizedbincount.c, na.rm = T)) * length(normalizedbincount.c)

# predict
enet.ma = bincounts.c %*% elnetbeta + elnetintercept
enet.c = enet.ma[1, 1]
wrsc.c = ff.pred(bincounts.c, B, mu, parameter.1, parameter.2)

# output
message(enet.c, '\t', wrsc.c, '\t', mean(c(enet.c, wrsc.c)))
temp.df = data.frame(file_name = BAM, enet = enet.c, wrsc = wrsc.c, avg = mean(c(enet.c, wrsc.c)))

output_name = paste0(basename(BAM), '.ff')
if (file.exists(output_name)) {
   i = 1
   while (TRUE) {
      output_name = paste0(basename(BAM), '_', i, '.ff')
      if (!file.exists(output_name)) {
         break
      }
      i = i + 1
   }
}

write_tsv(temp.df, output_name, col_names = F)



