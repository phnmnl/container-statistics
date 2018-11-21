#!/usr/bin/env Rscript

# Check if the limma and argparse are available, if not then install them
list.of.packages <- c("argparse", "limma","qvalue")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load argparse for parsing arguments
suppressWarnings(library("argparse"))

# set require arguments
parser <- ArgumentParser()

# set require arguments.
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")

parser$add_argument("-q", "--quietly", action="store_false",
                    dest="verbose", help="Print little output")

parser$add_argument("-in", "--dataMatrix_in", type="character",
                    help="Input data matrix containing peaks")

parser$add_argument("-s", "--sampleMetadata_in", type="character",
                    help="Input data matrix containing sample metadata")


parser$add_argument("-p", "--sampleVariable_in", type="character",
                    help="Input multiple data matrix containing variable information. They have to be separated by space or , or \\| or tab or |. The order MUST be identical to dataMatrix_in")


parser$add_argument("-out", "--sampleVariable_out", type="character",
                    help="Output data matrix containing peaks")


parser$add_argument("-c", "--covariates", type="character",
                    help="one or more of the column names (first row) of your sample metadata showing covariates you want to correct for. If more than one column selected, it should be separated by space or , or \\| or tab or |")


parser$add_argument("-g", "--group", type="character",
                    help="Treatment conditions to be preserved")


parser$add_argument("-f", "--contrasts", type="character",
                    help="One or more contrasts in the form of X-Y where levels of X and Y have to be in column specified using --group. If more than one contrasts selected, they should be separated by space or , or \\| or tab or |. Example, comparing groups X to Y and X to Z: X-Y,X-Z. No separating character in the contrast allowed!")

parser$add_argument("-b1", "--batch1", type="character",
                    help="one of the column names (first row) of your sample metadata showing batch variable (this will be factored)")

parser$add_argument("-b2", "--batch2", type="character",
                    help="one of the column names (first row) of your sample metadata showing batch variable (this will be factored)")

parser$add_argument("-r", "--regular", type="logical",default=FALSE,
                    help="If set, regular t-statistics will be performed")

# parse arguments
args <- parser$parse_args()


if(is.null(args$sampleVariable_out))
{
  errorMessage<-"No output file has been specified. You MUST specify the output file see the help (-h)!"
  write(errorMessage,stderr())
  stop(errorMessage,
       call. = FALSE)
}

if(is.null(args$sampleMetadata_in))
{
  errorMessage<-"No metadata input file has been specified. You MUST specify the input file see the help (-h)!"
  write(errorMessage,stderr())
  stop(errorMessage,
       call. = FALSE)
}

if(is.null(args$sampleMetadata_in))
{
  errorMessage<-"No metadata input file has been specified. You MUST specify the input file see the help (-h)!"
  write(errorMessage,stderr())
  stop(errorMessage,
       call. = FALSE)
}

if(is.null(args$sampleVariable_in))
{
  errorMessage<-"No variable input file has been specified. You MUST specify the input file see the help (-h)!"
  write(errorMessage,stderr())
  stop(errorMessage,
       call. = FALSE)
}

if(is.null(args$group))
{
  errorMessage<-"No group information has been specified. You MUST specify the group. see the help (-h)!"
  write(errorMessage,stderr())
  stop(errorMessage,
       call. = FALSE)
}




if(is.null(args$contrasts))
{
  errorMessage<-"No contrasts information has been specified. The script will assume one group test."
  write(errorMessage, stdout())
}



if ( args$verbose ) {
  write("Loading data matrix...\n", stdout())
}

# Load peak matrix

xMN <- t(as.matrix(read.table(args$dataMatrix_in,
                              check.names = FALSE,
                              header = TRUE,
                              row.names = 1,
                              sep = "\t",
                              comment.char = "")))


if ( args$verbose ) {
  write("Loading sample data...\n", stdout())
}

# Load sample metaData
samDF <- read.table(args$sampleMetadata_in,
                    check.names = FALSE,
                    header = TRUE,
                    row.names = 1,
                    sep = "\t",
                    comment.char = "")



if ( args$verbose ) {
  write("Loading variables ...\n", stdout())
}

# Load sample metaData
varDF <- read.table(args$sampleVariable_in,
                    check.names = FALSE,
                    header = TRUE,
                    row.names = 1,
                    sep = "\t",
                    comment.char = "")



# generate error message if row and column names are not identical
if(!identical(rownames(xMN), rownames(samDF)))
{
  errorMessage<-"Sample names (or number) in the data matrix (first row) and sample metadata (first column) are not identical; use the 'Check Format' module in the 'Quality Control' section "
  errorMessage<-paste(errorMessage,"files: ", dataMatrix_inFile,", ",sampleMetadata_inFile,"\n",sep="")
  write(errorMessage,stderr())
  stop(errorMessage,
       call. = FALSE)
}

if(!identical(colnames(xMN), rownames(varDF)))
{
  errorMessage<-"Variable names (or number) in the data matrix (first column) and sample variable (first row) are not identical; use the 'Check Format' module in the 'Quality Control' section "
  errorMessage<-paste(errorMessage,"files: ", dataMatrix_inFile,", ",sampleVariable_inFile,"\n",sep="")
  write(errorMessage,stderr())
  stop(errorMessage,
       call. = FALSE)
}


batch1<-NULL
# check if batch 1 is in meta data
if(!is.null(args$batch1) && args$batch1%in%colnames(samDF))
{
  batch1 <- matrix(samDF[, args$batch1], ncol = 1, dimnames = list(rownames(xMN), args$batch1))
  
}else if(!is.null(args$batch1)){
  
  errorMessage<-"Batch1 was not found in the column names (first row) of your sample metadata!"
  write(errorMessage,stderr())
  stop(errorMessage,
       call. = FALSE)
}

batch2<-NULL
# check if batch 2 is in meta data

if(!is.null(args$batch2) && args$batch2%in%colnames(samDF))
{
  batch2 <- matrix(samDF[, args$batch2], ncol = 1, dimnames = list(rownames(xMN), args$batch2))
}else if(!is.null(args$batch2)){
  
  
  errorMessage<-"Batch2 was not found in the column names (first row) of your sample metadata!"
  write(errorMessage,stderr())
  stop(errorMessage,
       call. = FALSE)
}


# for variable, first split it and then check all of them
covariates<-NULL

covariatesColumns<-ifelse(!is.null(args$covariates),
                          yes = sapply(strsplit(x = args$covariates,split = "\\;|,| |\\||\\t"),function(x){x}),
                          no = NA)

if(!is.null(args$covariates))
{
  covariatesColumns<-sapply(strsplit(x = args$covariates,split = "\\;|,| |\\||\\t"),function(x){x})
}else{covariates<-NULL}

if(!is.na(covariatesColumns) && all(covariatesColumns%in%colnames(samDF)))
{
  covariates <- samDF[, covariatesColumns]
}else if(!is.na(covariatesColumns)){
  
  notFoundColumns<-paste(covariatesColumns[!covariatesColumns%in%colnames(samDF)],sep=", ")
  errorMessage<-paste(notFoundColumns,"was/were not found in the column names (first row) of your sample metadata!")
  write(errorMessage,stderr())
  stop(errorMessage,
       call. = FALSE)
}else{
  covariates<-NULL
}


group<-NULL
design<-NULL
# check if batch 1 is in meta data
if(!is.null(args$group) && args$group%in%colnames(samDF))
{
  group <- matrix(samDF[, args$group], ncol = 1, dimnames = list(rownames(xMN), args$group))
  
  
}else if(!is.null(args$group)){
  
  errorMessage<-"Group was not found in the column names (first row) of your sample metadata!"
  write(errorMessage,stderr())
  stop(errorMessage,
       call. = FALSE)
}else{
  
  group<-rep(1,nrow(samDF))
}

design<-model.matrix(~0+group)
colnames(design)<-gsub(colnames(design),pattern = "group",replacement = "",fixed = T)
if ( args$verbose ) {
  write("Loading limme package ...\n", stdout())
}
suppressWarnings(library(limma))

X.batch<-NULL
if (!is.null(batch1)) {
  batch1 <- as.factor(batch1)
  contrasts(batch1) <- contr.sum(levels(batch1))
  batch1 <- model.matrix(~batch1)[, -1, drop = FALSE]
  colnames(batch1)<-args$batch1
}
if (!is.null(batch2)) {
  batch2 <- as.factor(batch2)
  contrasts(batch2) <- contr.sum(levels(batch2))
  batch2 <- model.matrix(~batch2)[, -1, drop = FALSE]
  colnames(batch2)<-args$batch1
}
if (!is.null(covariates)) 
  covariates <- as.matrix(covariates)
X.batch <- cbind(batch1, batch2, covariates)

if ( args$verbose ) {
  write("Creating design matrix ...\n", stdout())
}

design<-cbind(design, X.batch)



if ( args$verbose ) {
  write("Design matrix was created. Now creating contrast matrix ...\n", stdout())
}



design.pairs <-
  function(levels) {
    n <- length(levels)
    design <- matrix(0,n,choose(n,2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n,2)
    k <- 0
    for (i in 1:(n-1))
      for (j in (i+1):n) {
        k <- k+1
        design[i,k] <- 1
        design[j,k] <- -1
        colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
      }
    design
  }

contrasts<-NULL

if(!is.null(args$contrasts))
{
  contrasts<-sapply(strsplit(x = args$contrasts,split = "\\;|,| |\\||\\t"),function(x){x})
}

ContrastLevels<-NULL
designC<-NULL

if(!is.na(contrasts))
{
  ContrastLevels<-unique(sapply(strsplit(x = contrasts,split = "-",fixed=T),function(x){x}))
  
  if ( args$verbose ) {
    write("Checking if all the levels are in group...\n", stdout())
  }

  if(!all(ContrastLevels%in%group))
  {
    errorMessage<-"Not all the contrast levels are found in group column"
    write(errorMessage,stderr())
    stop(errorMessage,
         call. = FALSE)
  }
  designC<-(makeContrasts(contrasts=contrasts, levels=design)) 
  if ( args$verbose ) {
    write("Contrast matrix was created. Now preparing for statistics ...\n", stdout())
  }
  
}

if(!is.null(designC))
{
  fitF<-lmFit(t(xMN),design = design)
  fitS <- contrasts.fit(fit = fitF, contrasts = designC)
  
}else{
  
  fitS<-lmFit(t(xMN))
  
}

fitS <- eBayes(fitS)

convertLimmaToregularTtest<-function(fit)
{
  require(limma)
  ebfit <- eBayes(fit)
  ebfit$metabolites <- ebfit$genes
  ebfit$genes <- ebfit$rank <- ebfit$assign <- NULL
  ebfit$qr <- ebfit$qraux <- ebfit$pivot <- ebfit$tol <- NULL
  ebfit$cov.coefficients <- ebfit$pivot <- ebfit$lods <- NULL
  beta <- ebfit$coeff
  tstat <- sweep(data.matrix(ebfit$coef/ebfit$stdev.unscaled), 
                 1, ebfit$sigma, "/")
  tpval <- 2 * pt(-abs(tstat), df = ebfit$df.residual)
  ebfit$t <- tstat
  df <- ebfit$df.residual
  fstat <- classifyTestsF(ebfit, df = df, fstat.only = TRUE)
  Fstat <- as.vector(fstat)
  df1 <- attr(fstat, "df1")
  df2 <- attr(fstat, "df2")
  if (df2[1] > 1e+06) {
    Fpval <- pchisq(df1 * Fstat, df1, lower.tail = FALSE)
  }
  else {
    Fpval <- pf(Fstat, df1, df2, lower.tail = FALSE)
  }
  se <- beta/tstat
  ebfit$F <- Fstat
  ebfit$F.p.value <- Fpval
  ebfit$p.value <- tpval
  ebfit$df.total <- ebfit$s2.post <- ebfit$stdev.unscaled <- NULL
  ebfit$var.prior <- ebfit$proportion <- ebfit$s2.prior <- NULL
  ebfit$df.prior <- NULL
  ebfit$std.error <- se
  ebfit$t <- tstat
  ebfit$df <- df
  
  return(ebfit)
}
if(args$regular)
{
 write("Calculate regular t-statistics ...\n", stdout())
fitS<-convertLimmaToregularTtest(fitS)
}
fitSD<-data.frame(fitS)

if ( args$verbose ) {
  write("Calculate q-value ...\n", stdout())
}





fitSD$qvalue<-qvalue::qvalue(fitSD$F.p.value)$qvalue
varDF<-cbind(varDF,fitSD)


varDF<-cbind.data.frame(variableMetadata=rownames(varDF),varDF,stringsAsFactors = F)
if ( args$verbose ) {
  write("Writing output ...\n", stdout())
}

write.table(x = varDF,file = args$sampleVariable_out,
            row.names = F,quote = F,sep = "\t")


if ( args$verbose ) {
  write("Program finished!\n", stdout())
}


