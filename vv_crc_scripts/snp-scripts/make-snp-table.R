library(getopt)
library(trellis)


opts = getopt(matrix(c('normal','n',2,'character',
                       'tumor','t',2,'character',
                       'capture','c',2,'character',
                       'assembly','a',2,'character',
		       'output','o',2,'character'), ncol = 4, byrow = T))

sources = list(hg18 = 'svfilters.hg18', hg19 = 'svfilters.hg19', hg38 = 'svfilters.hg38')
package = sources[[opts$assembly]]
# data(snps, package = package)
data(dbsnp150_snps, package = package)

capture = read.table(opts$capture, header = FALSE, sep = '\t', stringsAsFactor = FALSE)
capture.ranges = GRanges(seqnames = capture[,1], ranges = IRanges(capture[,2], capture[,3]))

snps = svAF(normalBam = ops$normal,
            tumorBam = opts$tumor,
            genome = opts$assembly,
            n = length(dbsnp150_snps),
            position = dbsnp150_snps,
            region = capture.ranges)


snps = data.frame(snps)
out = snps[,c('seqnames','end','RefBase','AltBase','Normal.Mut.Count',
              'Normal.Coverage', 'Tumor.Mut.Count','Tumor.Coverage', 'Tumor.MAF')]
colnames(out)[1:2] = c('Chrom', 'Pos')
write.table(out, file  = opts$output, sep = '\t', quote = FALSE, row.names = FALSE)
