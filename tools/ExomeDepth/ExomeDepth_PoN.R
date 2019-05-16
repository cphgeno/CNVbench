library(optparse)
library(GenomicRanges)
library(ExomeDepth)
data(exons.hg19)
data(Conrad.hg19)

# Parse command line options
option_list <- list(make_option(c('-r','--ref'), action='store', type='character', default=NULL, help='Reference sequence (FASTA).'),
                    make_option(c('-b','--bed'), action='store', type='character', default=NULL, help='Regions file (BED).'),
                    make_option(c('-p','--pon'), action='store', type='character', default=NULL, help='File (list one per line) or folder with BAM files to be used as PoN.'),
                    make_option(c('-o','--out'), action='store', type='character', default=NULL, help='Output folder.'),
                    make_option(c('--debug'), action='store_true', type='logical', default=FALSE, help='Debug mode.')
)
opt <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE)
bam_files = opt$args
opt <- opt$options

if(opt$debug) {
  print(opt$ref)
  print(opt$bed)
  print(opt$pon)
  print(opt$bam_files)
}

# Parse BED file
regions <- read.table(opt$bed, header=FALSE, stringsAsFactors=FALSE)
colnames(regions) <- c("chromosome", "start", "end", "name")
# Get annotations
regions.GRanges <- GRanges(seqnames = regions$chromosome,
                           IRanges(start=regions$start,end=regions$end),
                           names = regions$name)

#######################
### Process samples ###
#######################
cat("Processing sample(s) ...", fill=TRUE)
# Get counts
bam_counts_gr <- GRanges(getBamCounts(bam.files=bam_files, bed.frame=regions, referenceFasta=opt$ref))
bam_names <- colnames(mcols(bam_counts_gr))[-1]
# Convert to data.frame
bam_counts_df <- data.frame(bam_counts_gr)
bam_counts_df$names <- names(bam_counts_gr)
# Fix column names
colnames(bam_counts_df) <- sub("[.]","-",sub("[.]","-",colnames(bam_counts_df)))
### prepare the main matrix of read count data
bam_counts <- as.matrix(subset(bam_counts_df, select=-c(seqnames, start, end, width, strand, GC, names)))
n_samples <- ncol(bam_counts)

###########################
### Process PoN samples ###
###########################
cat("Processing PoN ...", fill=TRUE)
# Get counts
if(dir.exists(opt$pon)) {
  ref_bam <- dir(opt$pon, "*.bam", full.names=TRUE)
} else if(file.exists(opt$pon)) {
  ref_bam <- read.table(opt$pon, header=FALSE, stringsAsFactors=FALSE)[1]
} else {
  stop("invalid PoN provided!")
}
ref_counts_gr <- GRanges(getBamCounts(bam.files=ref_bam, bed.frame=regions, referenceFasta=opt$ref))
# Convert to matrix
ref_counts <- as.matrix(subset(data.frame(ref_counts_gr), select=-c(seqnames, start, end, width, strand, GC)))
n_ref <- ncol(ref_counts)

### start looping over each sample
for (i in 1:n_samples){
  #### Create the aggregate reference set for this sample
  my.choice <- select.reference.set (test.counts =  bam_counts[,i],
                                     reference.counts = ref_counts,
                                     bin.length = bam_counts_df$width/1000,
                                     n.bins.reduced = 10000)
  
  my.reference.selected <- apply(X = ref_counts[, my.choice$reference.choice, drop = FALSE], MAR=1, FUN=sum)

  message('Now creating the ExomeDepth object')
  all.exons <- new('ExomeDepth',
                   test = bam_counts[,i],
                   reference = my.reference.selected,
                   formula ='cbind(test, reference) ~ 1')
  
  ################ Now call the CNVs
  all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        chromosome = bam_counts_df$seqnames,
                        start = bam_counts_df$start,
                        end = bam_counts_df$end,
                        name = bam_counts_df$names)
  
  ########################### Now annotate the ExomeDepth object
  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = Conrad.hg19.common.CNVs,
                             min.overlap = 0.5,
                             column.name ='Conrad.hg19')
  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = regions.GRanges,
                             min.overlap = 0.0001,
                             column.name ='regions')

  # Create output folder
  dir.create(opt$out, showWarnings = FALSE)
  # Write results to file  
  write.table(all.exons@CNV.calls,
              paste(opt$out, paste(sub("[.].*","",bam_names[i]), "tsv", sep="."), sep="/"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  # Plot results
#  plot(all.exons, count.threshold = 20, cex.lab = 0.8)
}
