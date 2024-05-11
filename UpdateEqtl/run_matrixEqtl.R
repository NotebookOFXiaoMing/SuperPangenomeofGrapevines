library(MatrixEQTL)

## 1  012 matrix
## 2  quant.norm tpm
## 3  covariates
## 4  snp pos
## 5  gene pos
args <- commandArgs(trailingOnly = TRUE)

SNP_file_name<- args[1]
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 20000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

expression_file_name<- args[2]

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

covariates_file_name<- args[3]

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$LoadFile(covariates_file_name);

snps_location_file_name<- args[4]
gene_location_file_name<- args[5]

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);



cisDist<-5000
pvOutputThreshold_cis<-0.1
pvOutputThreshold_tra<-0.00000000001

errorCovariance<-numeric()

useModel<-modelLINEAR

output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);


unlink(output_file_name_tra);
unlink(output_file_name_cis);

save(me,file = "me.Rdata")