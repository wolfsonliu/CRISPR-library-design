suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

source('parameters.r')

####################
con <- DBI::dbConnect(RSQLite::SQLite(), db.file)

spacers <- DBI::dbFetch(DBI::dbSendQuery(
    con,
    'SELECT * FROM gene_spacer'
))

## correct specific score

spacers$specific_score[is.na(spacers$specific_score)] <- 0.5

spacers$specific_score.corrected <- spacers$specific_score

spacers$specific_score.corrected[
            spacers$specific_score > 1
        ] <- 1 / ((1 / spacers$specific_score[spacers$specific_score > 1]) + 1)


####################
## distribution of features

## GC

p <- spacers %>%
    ggplot(aes(x=gc)) +
    geom_histogram(
        fill=RColorBrewer::brewer.pal(3, 'Set2')[2]
    ) +
    theme_classic()

ggsave('hist_gc.pdf', p)

p <- spacers %>%
    ggplot(aes(x=specific_score.corrected)) +
    geom_histogram(
        fill=RColorBrewer::brewer.pal(3, 'Set2')[2]
    ) +
    theme_classic()

ggsave('hist_specific_score.pdf', p)

p <- spacers %>%
    ggplot(aes(x=mismatch_score)) +
    geom_histogram(
        fill=RColorBrewer::brewer.pal(3, 'Set2')[2]
    ) +
    theme_classic()

ggsave('hist_mismatch_score.pdf', p)


####################

spacers <- spacers %>%
    group_by(gene_id) %>%
    arrange(exon_rank, specific_score)

spacers$punish <- (
    (
        as.numeric(spacers$exon_rank > 10) +
        as.numeric(spacers$exon_rank > 5) +
        as.numeric(spacers$exon_rank > 3)
    ) / 2
) +
    (5 - cut(spacers$specific_score.corrected, breaks=5, labels=FALSE)) +
    ((cut(abs(spacers$gc - 0.5), breaks=4, labels=FALSE) - 1) / 2)


p <- spacers %>%
    ggplot(aes(x=punish)) +
    geom_histogram(
        fill=RColorBrewer::brewer.pal(3, 'Set2')[2]
    ) +
    theme_classic()

ggsave('hist_punish.pdf', p)


genes <- unique(spacers$gene_id)

lib <- Reduce(
    rbind,
    spacers %>%
    group_by(gene_id) %>%
    group_map(
        ~ .x %>% transform(gene_id = .y) %>%
            arrange(punish) %>% head(3)
    )
)

gene.spacer.count <- lib %>%
    group_by(gene_id) %>% summarize(count=n())

lib$specific_score <- lib$specific_score.corrected

lib$specific_score.corrected <- NULL

write.table(
    lib[
        order(lib$gene_id),
        c('spacer', 'gene_id',
          'seqnames', 'start', 'end', 'strand',
          'pam', 'exon_rank', 'gc',
          'm0', 'm1', 'm2', 'm3',
          'specific_score', 'mismatch_score', 'punish')
    ],
    'library.txt',
    sep='\t', row.names=FALSE, quote=FALSE
)

####################
DBI::dbDisconnect(con)
####################
