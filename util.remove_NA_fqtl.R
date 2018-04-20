

source('Util.R')
library(readr)
library(dplyr)


clean.na <- function(chr) {

    tis.files <- list.files('result/fqtl_obs/chr' %&&% chr,
                            pattern = 'tis-factor.gz', full.names = TRUE)

    test.file <- function(ff) {
        temp <- read_tsv(ff, col_types = 'ciiddd') %>%
            na.omit()
        if(nrow(temp) < 1) {
            log.msg('Found NA %s', ff)
            return(ff)
        }
        return(NULL)
    }

    rm.files <- tis.files %>% sapply(test.file)
    rm.files <- do.call(c, rm.files) %>% .unlist()

    if(length(rm.files) > 0) {
        sapply(rm.files, gsub, pattern = '.tis-factor.gz', replacement = '.combined.gz') %>%
            sapply(FUN = file.remove)
        sapply(rm.files, gsub, pattern = '.tis-factor.gz', replacement = '.snp-factor.gz') %>%
            sapply(FUN = file.remove)
        sapply(rm.files, gsub, pattern = '.tis-factor.gz', replacement = '.snp-max.gz') %>%
            sapply(FUN = file.remove)
        sapply(rm.files, FUN = file.remove)
    }
}

sapply(seq(22, 1, -1), clean.na)
