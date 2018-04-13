
read.tis.sample.info <- function(fam.file = 'geno/chr1.fam', tis.file = 'data/tissues.txt') {

    require(dplyr)
    require(readr)    

    plink.fam <- read_delim(fam.file, delim = ' ', col_names = 'SUBJID', col_types = '_c____') %>%
        mutate(geno.pos = 1:n()) %>%
            mutate(SUBJID = gsub(SUBJID, pattern = 'GTEX[-]', replacement = ''))

    tis.tab <- read_tsv(tis.file, col_names = c('tis.dir', 'tis.name', 'tis.size'),
                        col_types = 'cci_') %>%
                            mutate(tis.idx = 1:n())

    sample.files <- 'data/' %&&% tis.tab$tis.dir %&&% '/samples.txt.gz'

    .read.samp.info <- function(...) {
        col.names <- c('SAMPID', 'SUBJID', 'data.pos', 'SMRIN', 'SMTSISCH',
                       'tis.name', 'SEX', 'AGE', 'HGHT', 'WGHT', 'BMI', 'SEX.std', 'AGE.std',
                       'HGHT.std', 'WGHT.std', 'BMI.std', 'ntis', 'tis.idx', 'tis.dir')
        ret <- suppressMessages(read_tsv(..., col_names = col.names) %>%
                                    mutate(expr.pos = 1:n()) %>%
                                        left_join(plink.fam))
        return(ret)
    }

    ret <- sample.files %>% lapply(.read.samp.info) %>%
        bind_rows()

    return(ret)
}
