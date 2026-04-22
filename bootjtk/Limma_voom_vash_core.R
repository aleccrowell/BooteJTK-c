#!/usr/bin/env Rscript
# Statistical core: f_changeNA + vooma (or voom) + vash shrinkage.
#
# Reads a pre-processed TSV produced by limma_preprocess.prepare_timeseries()
# (numeric column names, deduplicated row names) and writes a long-format TSV
# with columns: ID, Time, Mean, SD, SDpre, N.
#
# Usage:
#   Rscript Limma_voom_vash_core.R <in.tsv> <out.tsv> <period> [rnaseq]
library('limma')
library('vashr')

args        = commandArgs(trailingOnly=TRUE)
fn_in       = args[1]
fn_out      = args[2]
period      = as.numeric(args[3])
bool.rnaseq = length(args) > 3

df = read.csv(fn_in, header=TRUE, sep='\t', check.names=FALSE)
rownames(df) = df[, 1]
df = df[, -1]
tx = as.numeric(colnames(df))

t2  = sort(unique(tx %% period))
MAX = max(table(tx %% period))

f_changeNA <- function(x) {
  sd_val = median(apply(x, 1, function(y) sd(y, na.rm=TRUE)), na.rm=TRUE)
  m      = mean(apply(x, 2, function(z) mean(z, na.rm=TRUE)), na.rm=TRUE)
  sd_m   = sd(apply(x, 2, function(z) mean(z, na.rm=TRUE)), na.rm=TRUE)
  f_nareplace <- function(z) {
    if (is.na(mean(z, na.rm=TRUE))) {
      z <- rnorm(length(z), m, sd_m)
    } else if (sum(is.na(z)) > 0) {
      z[is.na(z)] <- rnorm(sum(is.na(z)), mean(z, na.rm=TRUE), sd_val)
    }
    return(z)
  }
  t(apply(x, 1, f_nareplace))
}

getmode <- function(v) {
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
f_isna <- function(x) sum(!is.na(x))

series.new  = NULL
times       = c()
rownames.id = rep(rownames(df), length(t2))

for (h in t2) {
  times = c(times, rep(h, nrow(df)))
}

for (h in t2) {
  cols = which(tx %% period == h)

  if (length(cols) == 1) {
    E_h = as.matrix(df[, cols, drop=FALSE])
    W_h = matrix(NaN, nrow=nrow(df), ncol=MAX)
    for (i in seq_len(MAX - 1)) {
      E_h = cbind(E_h, rep(NaN, nrow(df)))
    }
    ser.voom = list(E=E_h, weights=W_h)
  } else {
    ser = df[, cols]
    ser = f_changeNA(ser)
    if (bool.rnaseq) {
      ser.voom = voom(ser)
    } else {
      ser.voom = vooma(ser)
    }
    n_pad = MAX - ncol(ser.voom$E)
    for (i in seq_len(n_pad)) {
      ser.voom$weights = cbind(ser.voom$weights, rep(NaN, nrow(ser.voom$weights)))
      ser.voom$E       = cbind(ser.voom$E,       rep(NaN, nrow(ser.voom$E)))
    }
  }

  if (is.null(series.new)) {
    series.new = ser.voom
  } else {
    colnames(ser.voom$weights) = colnames(series.new$weights)
    series.new$weights = rbind(series.new$weights, ser.voom$weights)
    colnames(ser.voom$E) = colnames(series.new$E)
    series.new$E         = rbind(series.new$E, ser.voom$E)
  }
}

sds.pre = 1 / sqrt(series.new$weights[, 1])
sds.pre[is.na(sds.pre)] = runif(sum(is.na(sds.pre)), 1, length(sds.pre))
df.vash  = getmode(apply(series.new$E, 1, f_isna))
sds.vash = vash(sds.pre, df=df.vash)

fmean = function(x) mean(x, na.rm=TRUE)
f_n   = function(x) sum(!is.na(x))

result = data.frame(
  ID    = rownames.id,
  Time  = times,
  Mean  = apply(series.new$E, 1, fmean),
  SD    = sds.vash$sd.post,
  SDpre = sds.pre,
  N     = apply(series.new$E, 1, f_n),
  stringsAsFactors = FALSE
)

write.table(result, file=fn_out, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
