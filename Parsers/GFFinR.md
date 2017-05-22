# GFF Parser in R
* The GFF file format stands for "Gene Finding Format" or or "General Feature Format" and was invented at the Sanger Centre.
```r
getAttributeField <- function (x, field, attrsep = ";") {
     s = strsplit(x, split = attrsep, fixed = TRUE)
     sapply(s, function(atts) {
         a = strsplit(atts, split = "=", fixed = TRUE)
         m = match(field, sapply(a, "[", 1))
         if (!is.na(m)) {
             rv = a[[m]][2]
         }
         else {
             rv = as.character(NA)
         }
         return(rv)
     })
}
```
and here is quick parser
```r
gffRead <- function(gffFile, nrows = -1) {
     cat("Reading ", gffFile, ": ", sep="")
     gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
     header=FALSE, comment.char="#", nrows = nrows,
     colClasses=c("character", "character", "character", "integer",
"integer",
     "character", "character", "character", "character"))
     colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
     cat("found", nrow(gff), "rows with classes:",
         paste(sapply(gff, class), collapse=", "), "\n")
     stopifnot(!anyis.na(gff$start)), !anyis.na(gff$end)))
     return(gff)
}
```
Now you can do stuff like
```r
gff <- gffRead(gfffile)
gff$Name <- getAttributeField(gff$attributes, "Name")
gff$ID <- getAttributeField(gff$attributes, "ID")
```
gfffile is just an object holding the file name.
