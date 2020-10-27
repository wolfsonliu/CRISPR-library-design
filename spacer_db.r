if (!requireNamespace('BiocManager', quietly = TRUE)) {
    install.packages('BiocManager')
}
if (!requireNamespace('Biostrings', quietly = TRUE)) {
    BiocManager::install('Biostrings')
}
if (!requireNamespace('GenomicFeatures', quietly = TRUE)) {
    BiocManager::install('GenomicFeatures')
}
if (!requireNamespace('GenomicRanges', quietly = TRUE)) {
    BiocManager::install('GenomicRanges')
}
if (!requireNamespace('RSQLite', quietly = TRUE)) {
    install.packages(c('DBI', 'RSQLite'))
}

####################

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
options(stringsAsFactors=FALSE)

runBowtie <- function(reference,
                      query,
                      output,
                      threads=1,
                      bowtie='bowtie') {
    ## -a, -n <core mismatches +1>, -l <core length>, -e <total mismatches  30 + 30> and -y.
    ## from CCtop
    system(
        paste(bowtie,'-a', '-n', '3', '-l', 20, '-e', '120', '-y',
              '-p', threads, '-f', reference, query, output)
    )
}

getSequence <- function(gr, chromosomes)
{                               # from motifRG package
    all.seq <- DNAStringSet(rep("", length(gr)))
    for(chr in (sort(unique(seqnames(gr))))){
        select <- as.vector(seqnames(gr)==chr)
        gr.select <- gr[select]
        chr.seq <- chromosomes[[chr]]
        seq <- DNAStringSet(Views(chr.seq, start=start(gr.select), end=end(gr.select)))
        neg <- as.vector(strand(gr.select) == "-" )
        seq[neg] <- reverseComplement(DNAStringSet(seq[neg]))
        all.seq[select] <- seq
    }
    all.seq
}

####################

source('parameters.r')

pam.pattern <- list(
    DNAString(pam),
    reverseComplement(DNAString(pam))
)

spacer.pattern <- DNAStringSet(
    c(
        xscat(
            paste(
                rep('N', times=spacer.length),
                collapse=''
            ),
            pam.pattern[[1]]
        ),
        xscat(
            pam.pattern[[2]],
            paste(
                rep('N', times=spacer.length),
                collapse=''
            )
        )
    )
)

####################
## reference data
transcript.anno <- makeTxDbFromGFF(
    gff.file,
    format='gff3',
    organism = species.name,
    taxonomyId = species.id,
    metadata=data.frame(
        name=c(
            'genome-build',
            'genome-build-accession',
            'resource url',
            'species url'
        ),
        value=c(
            genome.build,
            genome.build.accession,
            genome.build.url,
            species.url
        )
    )
)

saveDb(transcript.anno, TxDb.file)

transcript.anno <- loadDb(TxDb.file)

chromosomes <- readDNAStringSet(fa.file, format='fasta')

####################
## find all site with pam on the genome
con <- dbConnect(RSQLite::SQLite(), db.file)

dbCreateTable(
    con, 'spacer_pam_site',
    fields=c(
        'seqnames'='TEXT',
        'start'='INTEGER',
        'end'='INTEGER',
        'width'='INTEGER',
        'strand'='TEXT',
        'spacer'='TEXT',
        'pam'='TEXT',
        'gene_id'='TEXT',
        'exon_rank'='INTEGER'
    )
)

for (chrname in names(chromosomes)) {
    ## forward strand
    s.fwd <- matchPattern(
        spacer.pattern[[1]],
        chromosomes[[chrname]],
        fixed=FALSE
    )
    s.fwd.gr <- GRanges(
        seqnames=chrname,
        ranges=IRanges(
            start=start(s.fwd),
            end=end(s.fwd)
        ),
        strand='+',
        spacer=substr(s.fwd, 1, spacer.length),
        pam=substr(s.fwd, spacer.length + 1, spacer.length + pam.length),
        gene_id='',
        exon_rank=-1
    )
    print('Granges')
    hits <- findOverlaps(
        s.fwd.gr, exons(transcript.anno),
        select='first', ignore.strand=TRUE
    )
    print('find overlaps')
    mcols(s.fwd.gr[seq(length(hits))[!is.na(hits)]])$gene_id <- unlist(exons(transcript.anno, columns='gene_id')[na.omit(hits)]$gene_id)
    mcols(s.fwd.gr[seq(length(hits))[!is.na(hits)]])$exon_rank <- unlist(exons(transcript.anno, columns='exon_rank')[na.omit(hits)]$exon_rank)
    dbAppendTable(
        con, 'spacer_pam_site',
        as.data.frame(s.fwd.gr) %>%
        filter(pam %in% c('AGG', 'TGG', 'CGG', 'GGG')) %>%
        filter(!grepl('T{4}', spacer)) %>%
        filter(!grepl('A{10}', spacer)) %>%
        filter(!grepl('C{10}', spacer)) %>%
        filter(!grepl('G{10}', spacer))
    )
    print('append table')
    rm(s.fwd, s.fwd.gr)
    ## reverse strand
    s.rev <- matchPattern(
        spacer.pattern[[2]],
        chromosomes[[chrname]],
        fixed=FALSE
    )
    s.rev.gr <- GRanges(
        seqnames=chrname,
        ranges=IRanges(
            start=start(s.rev),
            end=end(s.rev)
        ),
        strand='-',
        spacer=substr(
            reverseComplement(s.rev),
            1, spacer.length
        ),
        pam=substr(
            reverseComplement(s.rev),
            spacer.length + 1,
            spacer.length + pam.length
        ),
        gene_id='',
        exon_rank=-1
    )
    hits <- findOverlaps(
        s.rev.gr, exons(transcript.anno),
        select='first', ignore.strand=TRUE
    )
    mcols(s.rev.gr[seq(length(hits))[!is.na(hits)]])$gene_id <- unlist(exons(transcript.anno, columns='gene_id')[na.omit(hits)]$gene_id)
    mcols(s.rev.gr[seq(length(hits))[!is.na(hits)]])$exon_rank <- unlist(exons(transcript.anno, columns='exon_rank')[na.omit(hits)]$exon_rank)

    ## save spacer sites to database table: spacer_pam_site
    dbAppendTable(
        con, 'spacer_pam_site',
        as.data.frame(s.rev.gr) %>%
        filter(pam %in% c('AGG', 'TGG', 'CGG', 'GGG')) %>%
        filter(!grepl('T{4}', spacer)) %>%
        filter(!grepl('A{10}', spacer)) %>%
        filter(!grepl('C{10}', spacer)) %>%
        filter(!grepl('G{10}', spacer))
    )
    rm(s.rev, s.rev.gr)
}

## get all unique spacer sequence
dbExecute(
    con,
    'CREATE TABLE IF NOT EXISTS spacer (
         spacer TEXT PRIMARY KEY,
         count  INTEGER,
         UNIQUE(spacer) ON CONFLICT IGNORE
    )'
)

dbExecute(
    con,
    'INSERT OR IGNORE INTO spacer(spacer, count)
     SELECT spacer, count(*)
     FROM spacer_pam_site
     GROUP BY spacer'
)

dbDisconnect(con)

## calculate spacer off target information

n.split <- 1000

## nbase <- ceiling(log(n.split, 4))

nbase <- 6

start.bases <- unlist(apply(
    expand.grid(
        rep(list(c('A', 'T', 'C', 'G')), times=nbase),
        stringsAsFactors=FALSE
    ),
    MARGIN=1,
    paste,
    collapse=''
))

con <- dbConnect(RSQLite::SQLite(), db.file)
dbCreateTable(
    con, 'spacer_off_target_raw',
    fields=c(
        'spacer'='TEXT',
        'mismatch'='INTEGER',
        'count'='INTEGER',
        'cfd'='REAL',
        'mismatch_score'='REAL'
    )
)
dbDisconnect(con)

con <- dbConnect(RSQLite::SQLite(), db.file)

load('mismatch.activate.Rdata')

if (parallel.threads == 1) {
    ## use non parallel strategy
    for (start.base in start.bases) {
        tmplabel <- paste(sample(c(letters, as.character(1:0)), 5), collapse='')
        tmpfa <- file.path(
            tmpdir,
            paste(tmplabel, 'fa', sep='.')
        )
        while (file.exists(tmpfa)) {
            tmplabel <- paste(sample(c(letters, as.character(1:0)), 5), collapse='')
            tmpfa <- file.path(
                tmpdir,
                paste(tmplabel, 'fa', sep='.')
            )
        }
        tmpoff <- file.path(
            tmpdir,
            paste(tmplabel, 'off', sep='.')
        )

        res <- dbSendQuery(
            con, 'SELECT spacer FROM spacer WHERE substr(spacer, 1, ?) = ? AND count = 1',
            params=list(
                nchar(start.base),
                start.base
            )
        )
        spacer <- dbFetch(res)
        if (dim(spacer)[1] > 0) {
            spacer <- DNAStringSet(spacer$spacer)
            names(spacer) <- as.character(spacer)
            writeXStringSet(spacer, tmpfa, format='fasta')

            runBowtie(
                index.prefix,
                tmpfa, tmpoff,
                threads=bowtie.threads
            )

            dat <- read.table(
                tmpoff, header=FALSE, sep='\t',
                stringsAsFactors=FALSE
            )
            file.remove(tmpfa)
            file.remove(tmpoff)
            gr <- GRanges(
                seqnames=dat$V3,
                ranges=IRanges(
                    start=ifelse(
                        dat$V2 == '+',
                        dat$V4 + spacer.length,
                        dat$V4 - pam.length
                    ),
                    end=ifelse(
                        dat$V2 == '+',
                        dat$V4 + spacer.length + pam.length -1,
                        dat$V4 - 1
                    )
                ),
                strand=dat$V2
            )

            mcols(gr)$pam <- getSequence(gr, chromosomes)
            selects <- vcountPattern(
                pam.pattern[[1]], gr$pam, fixed=FALSE
            ) > 0
            dat <- dat[selects,]
            dat$mismatch <- unlist(
                lapply(strsplit(dat$V8, ','), length)
            )
            ## CFD value
            dat$cfd <- unlist(
                lapply(
                    strsplit(dat$V8, ','),
                    function(a) {
                        if (length(a) == 0) {
                            return(1)
                        }
                        else {
                            return(Reduce(`*`, mismatch.activate[a]))
                        }
                    }
                )
            )
            ## CCtop mismatch_score
            dat$mismatch_score <- unlist(
                lapply(
                    strsplit(dat$V8, ','),
                    function(a) {
                        if (length(a) == 0) {
                            return(0)
                        }
                        else {
                            return(
                                sum(1.2 ^ as.numeric(gsub(':.*', '', a)))
                            )
                        }
                    }
                )
            )

            x <- dat %>%
                group_by(V1, mismatch, cfd, mismatch_score) %>%
                summarize(count=n()) %>%
                transform(
                    spacer=trimws(V1)
                ) %>%
                select(spacer, mismatch, count, cfd, mismatch_score)
            dbAppendTable(con, 'spacer_off_target_raw', x)
        }
    }

    dbDisconnect(con)
} else if (parallel.threads > 1) {
    cl <- makeCluster(parallel.threads)
    clusterExport(cl, "mismatch.activate")
    clusterExport(cl, "chromosomes")
    clusterExport(cl, "fa.file")
    clusterExport(cl, "db.file")
    clusterExport(cl, "index.prefix")
    clusterExport(cl, "spacer.length")
    clusterExport(cl, "pam.pattern")
    clusterExport(cl, "pam.length")
    clusterExport(cl, "tmpdir")
    clusterExport(cl, "bowtie.threads")
    clusterExport(cl, "runBowtie")
    clusterExport(cl, "getSequence")
    clusterEvalQ(cl, suppressPackageStartupMessages(library(Biostrings)))
    clusterEvalQ(cl, suppressPackageStartupMessages(library(dplyr)))
    clusterEvalQ(cl, suppressPackageStartupMessages(library(GenomicRanges)))

    tmpfiles <- parLapplyLB(
        cl,
        start.bases,
        function(start.base) {
            tmpfa <- file.path(
                tmpdir,
                paste(start.base, 'fa', sep='.')
            )
            tmpoff <- file.path(
                tmpdir,
                paste(start.base, 'off', sep='.')
            )
            con <- DBI::dbConnect(RSQLite::SQLite(), db.file)
            res <- DBI::dbSendQuery(
                            con,
                            'SELECT spacer FROM spacer WHERE substr(spacer, 1, ?) = ? AND count = 1',
                            params=list(
                                nchar(start.base),
                                start.base
                            )
                        )
            spacer <- DBI::dbFetch(res)
            DBI::dbDisconnect(con)
            if (dim(spacer)[1] > 0) {
                spacer <- DNAStringSet(spacer$spacer)
                names(spacer) <- as.character(spacer)
                writeXStringSet(spacer, tmpfa, format='fasta')

                runBowtie(
                    index.prefix,
                    tmpfa, tmpoff,
                    threads=bowtie.threads
                )

                dat <- read.table(
                    tmpoff, header=FALSE, sep='\t',
                    stringsAsFactors=FALSE
                )
                file.remove(tmpfa)
                file.remove(tmpoff)
                gr <- GRanges(
                    seqnames=dat$V3,
                    ranges=IRanges(
                        start=ifelse(
                            dat$V2 == '+',
                            dat$V4 + spacer.length,
                            dat$V4 - pam.length
                        ),
                        end=ifelse(
                            dat$V2 == '+',
                            dat$V4 + spacer.length + pam.length -1,
                            dat$V4 - 1
                        )
                    ),
                    strand=dat$V2
                )

                mcols(gr)$pam <- getSequence(gr, chromosomes)
                selects <- vcountPattern(
                    pam.pattern[[1]], gr$pam, fixed=FALSE
                ) > 0
                dat <- dat[selects,]
                dat$mismatch <- unlist(
                    lapply(strsplit(dat$V8, ','), length)
                )
                dat$cfd <- unlist(
                    lapply(
                        strsplit(dat$V8, ','),
                        function(a) {
                            if (length(a) == 0) {
                                return(1)
                            }
                            else {
                                return(Reduce(`*`, mismatch.activate[a]))
                            }
                        }
                    )
                )
                dat$mismatch_score <- unlist(
                    lapply(
                        strsplit(dat$V8, ','),
                        function(a) {
                            if (length(a) == 0) {
                                return(0)
                            }
                            else {
                                return(
                                    sum(1.2 ^ as.numeric(gsub(':.*', '', a)))
                                )
                            }
                        }
                    )
                )

                x <- dat %>%
                    group_by(V1, mismatch, cfd, mismatch_score) %>%
                    summarize(count=n()) %>%
                    transform(
                        spacer=trimws(V1)
                    ) %>%
                    select(spacer, mismatch, count, cfd, mismatch_score)

                tmpout <- file.path(
                    tmpdir, paste(start.base, 'spacer_off_target_raw', 'txt', sep='.')
                )
                write.table(
                    x,
                    tmpout,
                    row.names=FALSE, sep='\t', quote=FALSE
                )
                return(tmpout)
            }
        }
    )

    stopCluster(cl)

    con <- dbConnect(RSQLite::SQLite(), db.file)
    tmpfiles <- list.files(
        tmpdir, '.spacer_off_target_raw.txt'
    )
    for (x in tmpfiles) {
        xdata <- read.table(
            file.path(tmpdir, x),
            header=TRUE, sep='\t', stringsAsFactor=FALSE
        )

        dbAppendTable(con, 'spacer_off_target_raw', xdata)
    }
    dbDisconnect(con)
} else {
    stop('parallel.threads should be larger or equal to 1')
}

####################

con <- dbConnect(RSQLite::SQLite(), db.file)

dbExecute(
    con,
    'CREATE TABLE IF NOT EXISTS spacer_mismatch (
        spacer   TEXT,
        mismatch INTEGER,
        count    INTEGER,
        PRIMARY KEY (spacer, mismatch)
    );'
)

## add default counts to avoid the failure capture by bowtie
dbExecute(
    con,
    'INSERT INTO spacer_mismatch(spacer, mismatch, count)
     SELECT spacer, 0, count
     FROM spacer
     WHERE count = 1;'
)

dbExecute(
    con,
    'INSERT OR REPLACE INTO spacer_mismatch(spacer, mismatch, count)
     SELECT spacer, mismatch, sum(count)
     FROM spacer_off_target_raw
     GROUP BY spacer, mismatch;'
)

## table: spacer specific score
dbCreateTable(
    con, 'spacer_specific_score',
    fields=c(
        'spacer'='TEXT',
        'specific_score'='REAL'
    )
)

dbExecute(
    con,
    'INSERT INTO spacer_specific_score
         (spacer, specific_score)
     SELECT spacer, 1 / sum(count * cfd)
     FROM spacer_off_target_raw
     GROUP BY spacer;'
)

## table: spacer mismatch score
dbCreateTable(
    con, 'spacer_mismatch_score',
    fields=c(
        'spacer'='TEXT',
        'mismatch_score'='REAL'
    )
)

dbExecute(
    con,
    'INSERT INTO spacer_mismatch_score
         (spacer, mismatch_score)
     SELECT spacer,
     sum(mismatch_score * mismatch)
     FROM spacer_off_target_raw
     GROUP BY spacer;'
)

## table: spacer gc
dbCreateTable(
    con, 'spacer_gc',
    fields=c(
        'spacer'='TEXT',
        'gc'='REAL'
    )
)

spacer.gc <- dbFetch(dbSendQuery(
    con,
    'SELECT DISTINCT spacer FROM spacer_mismatch;'
))

spacer.gc$gc <- letterFrequency(
    DNAStringSet(spacer.gc$spacer),
    letters='GC', as.prob=TRUE
)[,1]

dbAppendTable(
    con, 'spacer_gc',
    spacer.gc[c('spacer', 'gc')]
)

dbDisconnect(con)

####################
## select gene targeting sgRNAs

con <- dbConnect(RSQLite::SQLite(), db.file)

spacer.off.target <- dbFetch(dbSendQuery(
    con, 'SELECT * FROM spacer_mismatch;'
)) %>%
    spread(key=mismatch, value=count, fill=0)

colnames(spacer.off.target) <- c('spacer', 'm0', 'm1', 'm2', 'm3')

spacer.gc <- dbFetch(dbSendQuery(
    con, 'SELECT * FROM spacer_gc;'
))

spacer.specific_score <- dbFetch(dbSendQuery(
    con, 'SELECT * FROM spacer_specific_score;'
))

spacer.mismatch_score <- dbFetch(dbSendQuery(
    con, 'SELECT * FROM spacer_mismatch_score;'
))

gene.targeting.spacer <- dbFetch(dbSendQuery(
    con,
    'SELECT * FROM spacer_pam_site
    WHERE length(gene_id) > 0;'
))

selected.spacer <- inner_join(
    spacer.off.target %>%
    filter(m0 <= 1) %>%
    filter(m1 == 0),
    gene.targeting.spacer,
    by='spacer'
)

selected.spacer <- selected.spacer[
    !(duplicated(selected.spacer$spacer) | duplicated(selected.spacer$spacer, fromLast = T)),
    ]

selected.spacer <- left_join(
    selected.spacer,
    spacer.gc,
    by='spacer'
)

selected.spacer <- left_join(
    selected.spacer,
    spacer.specific_score,
    by='spacer'
)

selected.spacer <- left_join(
    selected.spacer,
    spacer.mismatch_score,
    by='spacer'
)

## spacers targeting genes
dbExecute(
    con,
    'CREATE TABLE IF NOT EXISTS gene_spacer (
         spacer         TEXT PRIMARY KEY,
         seqnames       TEXT,
         start          INTEGER,
         end            INTEGER,
         strand         TEXT,
         pam            TEXT,
         gene_id        TEXT,
         exon_rank      INTEGER,
         gc             REAL,
         m0             INTEGER,
         m1             INTEGER,
         m2             INTEGER,
         m3             INTEGER,
         specific_score REAL,
         mismatch_score REAL,
         UNIQUE(spacer) ON CONFLICT IGNORE
    );'
)

dbAppendTable(
    con, 'gene_spacer',
    selected.spacer[
        c('spacer',
          'seqnames', 'start', 'end', 'strand',
          'pam', 'gene_id', 'exon_rank', 'gc',
          'm0', 'm1', 'm2', 'm3',
          'specific_score', 'mismatch_score')
    ]
)

dbDisconnect(con)

####################
