#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(plink2R)
})

option_list = list(
  make_option("--sumstats",    type="character", help="Summary‐stats file (SNP as CHR:BP, A1, A2, Z)"),
  make_option("--weights",     type="character", help="weights_list.txt"),
  make_option("--weights_dir", type="character", help="dir containing .wgt.RDat files"),
  make_option("--ref_ld_chr",  type="character", help="prefix to per-chr PLINK files (e.g. yeast1011_chr)"),
  make_option("--chr",         type="integer",   help="Chromosome number (1–16)"),
  make_option("--out",         type="character", help="Output prefix")
)
opt = parse_args(OptionParser(option_list=option_list))

#── 1) Load & parse summary stats ─────────────────────────────────────────────
ss = fread(opt$sumstats, sep="\t", header=TRUE, data.table=FALSE)
parts = tstrsplit(ss$SNP, ":", fixed=TRUE)
ss$CHRnum = as.integer(parts[[1]])
ss$BP     = as.integer(parts[[2]])
if (any(is.na(ss$CHRnum)|is.na(ss$BP)))
  stop("Cannot parse CHR:BP from summary‐stats")
ss = ss[ss$CHRnum == opt$chr, , drop=FALSE]
if (nrow(ss)==0) stop("No summary rows for chr ", opt$chr)

#── 2) Load & filter weights_list ────────────────────────────────────────────
wgtlist = read.table(opt$weights, header=TRUE, sep="\t", as.is=TRUE)
rom = as.character(as.roman(opt$chr))
wgtlist = wgtlist[ as.character(wgtlist$CHR) == rom, , drop=FALSE ]
if (nrow(wgtlist)==0) stop("No weight entries for chr ", opt$chr)
N = nrow(wgtlist)

#── 3) Load LD reference ─────────────────────────────────────────────────────
bed_pref = paste0(opt$ref_ld_chr, opt$chr)
genos = read_plink(bed_pref, impute="avg")
bim_bp = genos$bim[,4]

#── 4) Match summary→LD by BP ────────────────────────────────────────────────
m = match(bim_bp, ss$BP)
if (all(is.na(m))) stop("No positions match between sumstats and LD ref")
keep = !is.na(m)
genos$bim = genos$bim[keep, ]
genos$bed = genos$bed[,keep]
ss_sub   = ss[m[keep], ]
ss_sub$SNP = paste0(opt$chr, ":", ss_sub$BP)  # new SNP names

#── 5) Prepare output table ──────────────────────────────────────────────────
out = data.frame(
  FILE=character(N), ID=character(N), CHR=integer(N),
  P0=integer(N), P1=integer(N),
  NSNP=integer(N), NWGT=integer(N), MODEL=character(N),
  MODELCV.R2=character(N), MODELCV.PV=character(N),
  TWAS.Z=numeric(N), TWAS.P=numeric(N),
  stringsAsFactors=FALSE
)

#── 6) allele‐QC helper ──────────────────────────────────────────────────────
allele.qc = function(a1,a2,ref1,ref2) {
  a1=toupper(a1); a2=toupper(a2)
  ref1=toupper(ref1); ref2=toupper(ref2)
  comp = function(x) c(A="T",T="A",G="C",C="G")[x]
  flip1 = comp(ref1); flip2 = comp(ref2)
  keep = !((a1=="A"&a2=="T")|(a1=="T"&a2=="A")|
           (a1=="C"&a2=="G")|(a1=="G"&a2=="C"))
  keep = keep & a1 %in% c("A","T","G","C") & a2 %in% c("A","T","G","C")
  flip = (a1==ref2 & a2==ref1)|(a1==flip2 & a2==flip1)
  list(keep=keep, flip=flip)
}

#── 7) Main loop over weights files ─────────────────────────────────────────
for (i in seq_len(N)) {
  wf = wgtlist$WGT[i]
  load(file.path(opt$weights_dir, wf))   # loads snps, wgt.matrix, cv.performance

  # 7a) match snps→LD
  mm   = match(snps[,2], genos$bim[,2])
  ok   = !is.na(mm)
  snps      = snps[ok, ,drop=FALSE]
  wmat      = wgt.matrix[ok,,drop=FALSE]
  cur.bim   = genos$bim[mm[ok], ]
  cur.genos = scale(genos$bed[, mm[ok]])

  # 7b) flip alleles in weights if needed
  qc = allele.qc(snps[,5], snps[,6], cur.bim[,5], cur.bim[,6])
  qc$flip[is.na(qc$flip)] = FALSE
  wmat[qc$flip, ] = -wmat[qc$flip, ]

  # 7c) grab GWAS Z by BP
  m2 = match(cur.bim[,4], ss_sub$BP)
  Z  = ss_sub$Z[m2]
  Z[!qc$keep] = NA

  # 7d) best model from cv.performance
  best = which.min(cv.performance["pval", ,drop=TRUE])

  # 7e) TWAS Z & P
  twasz = sum(wmat[,best] * Z, na.rm=TRUE)
  twasv = as.numeric(wmat[,best] %*% cov(cur.genos) %*% wmat[,best])
  twasp = 2 * pnorm(-abs(twasz))

  # 7f) fill table row i
  out$FILE[i]       = wf
  out$ID[i]         = wgtlist$ID[i]
  out$CHR[i]        = opt$chr
  out$P0[i]         = wgtlist$P0[i]
  out$P1[i]         = wgtlist$P1[i]
  out$NSNP[i]       = nrow(cur.bim)
  out$NWGT[i]       = sum(wmat[,best]!=0)
  out$MODEL[i]      = colnames(wmat)[best]
  out$MODELCV.R2[i] = sprintf("%.3g", cv.performance["rsq",best])
  out$MODELCV.PV[i] = sprintf("%.3g", cv.performance["pval",best])
  out$TWAS.Z[i]     = twasz
  out$TWAS.P[i]     = twasp
}

#── 8) write out ─────────────────────────────────────────────────────────────
write.table(out, file=paste0(opt$out,".dat"),
            sep="\t", quote=FALSE, row.names=FALSE)

