#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)
library(broom)
library(jsonlite)

# Make sure we get reproducible results
set.seed(1234)

# We are doing three things in this script, seperated into three functions:
# 1. Simulate base covariates and a continuous phenotype <-  func simulate_base_covariates_and_phenotype
# 2. Simulate effects of common variants on phenotype <- func simulate_genetic_effects
# 3. Simulate rare WES variants <- func simulate_wes_data
simulate_base_covariates_and_phenotype <- function() {

  # First, we need to generate a set of base covariates
  # Input the PCs from our simulated data from simuPOP
  covars <- fread("sim_data/sim_chromosome.eigenvec", colClasses = c("IID"="character","FID"="character"))

  # Need to make a facsimile of the base_covariates.tsv file that is provided to all runassociationtesting runs:
  # Generate the following random covariates for all 10k participants; with column name in the file:
  # sex; 22001-0.0
  # age; 21003-0.0 (need to make sure to simulate a quadratic effect)
  # array batch; 22000-0.0
  # PCs; 22009-0.[1-40] (already done)
  # wes batch; wes.batch
  # genetics qc status; genetics_qc_pass
  covars[,eid:=IID]
  covars[,sex:=sample(c(0,1), size=nrow(covars), replace=T, prob=c(0.55, 0.45))]  # sex is slightly unbalanced in UKB
  covars[,age:=sample(c(40:70), size=nrow(covars), replace=T)]
  covars[,age_squared:=age ^ 2]
  covars[,array:=sample(c(-11:-1,1:95), size=nrow(covars), replace=T)]
  covars[,array:=if_else(array < 0, paste0('bileve',abs(array)), paste0('axiom',array))]
  covars[,wes.batch:=sample(c("200k", "450k", "470k"), size=nrow(covars), replace=T, prob=c(0.42, 0.53, 0.05))]
  covars[,genetics_qc_pass:=1]  # Have already tested this functionality, no need to do again here but still must be in the phenofile

  # And add a quantitative and categorical phenotype to test additional covariate functionality:
  covars[,quant_covar_1:=rnorm(nrow(covars), mean=2, sd=1)]
  covars[,cat_covar_2:=sample(c('bender', 'flexo', 'crushinator'),size=nrow(covars), replace=T, prob=c(0.6,0.3,0.1))]

  # Now we are starting to build a phenotype. Start with the PCs and generate a beta value for our phenotype
  # We want the 1st 10 PCs to explain ~10% of the total variance of the phenotype
  get.beta <- function(covar, var_desired) {
    # This runs the formula to scale individual PC beta based on desired variance
    # To test this and see the variance, do:
    # var(calc_beta * sim[,PC1])

    act_var <- var(covars[,get(covar)])
    calc_beta <- 10^((log10(var_desired) - log10(act_var)) / 2)
    return(calc_beta)

  }

  # Generate a data.table of PC
  sim_var_cont <- data.table("covar" = 1:10)
  sim_var_cont[,var_expl:=0.025*covar^-2]  # Basically, I want a decreasing amount of variance explained as PCs go on
  sim_var_cont[,covar:=paste0("PC",covar)]
  sim_var_cont <- rbind(sim_var_cont, data.table("covar"=c("sex","age","quant_covar_1","age_squared"), "var_expl"=c(0.02, 0.025, 0.04, 0.02)))
  sim_var_cont[,beta:=get.beta(covar,var_expl),by=seq_len(nrow(sim_var_cont))]
  sim_var_cont[,beta:=beta * c(sample(c(-1,1),12, replace = T), 1, 1)]  # Randomly simulate pos / neg effects for everything other than age

  # var.epsilon is just the cumulative variance and covariance of all the covariates in the model
  var.epsilon <- 0
  for (x in seq_len(nrow(sim_var_cont))) {
    x_beta <- sim_var_cont[x, beta]
    x_covar <- sim_var_cont[x, covar]
    var.epsilon <- var.epsilon + var(x_beta*covars[,get(x_covar)])
    for (y in seq_len(nrow(sim_var_cont))) {
      y_beta <- sim_var_cont[y, beta]
      y_covar <- sim_var_cont[y, covar]
      if (y <= x) { # Don't add cov(X,X) (which is just variance) or add cov(PCX,PCY) twice (as cov(PCY,PCX) == cov(PCX,PCY))
        next
      } else {
        var.epsilon <- var.epsilon + (2*cov(x_beta*covars[,get(x_covar)], y_beta*covars[,get(y_covar)]))
      }
    }
  }

  # Simulate categorical dummy variables the best we can and make sure they are included in our overall variance
  wes_betas <- data.table('wes.batch'=covars[,unique(wes.batch)], 'wes.betas'=rnorm(3, mean = 0, sd = 0.2))
  array_betas <- data.table('array'=covars[,unique(array)], 'array.betas'=rnorm(length(covars[,unique(array)]), mean = 0, sd = 0.05))
  cat_covar_betas <- data.table('cat_covar_2'=covars[,unique(cat_covar_2)], 'cat_covar_2.betas'=rnorm(3, mean=0,sd = 0.2))

  covars <- merge(covars, wes_betas, by='wes.batch')
  covars <- merge(covars, array_betas, by='array')
  covars <- merge(covars, cat_covar_betas, by='cat_covar_2')
  var.epsilon <- var.epsilon + var(covars[,wes.betas]) + var(covars[,array.betas]) + var(covars[,cat_covar_2.betas])
  var.epsilon <- var.epsilon + (2*cov(covars[,wes.betas], covars[,array.betas])) + (2*cov(covars[,wes.betas], covars[,cat_covar_2.betas])) + (2*cov(covars[,cat_covar_2.betas], covars[,array.betas]))

  for (y in seq_len(nrow(sim_var_cont))) {

    y_beta <- sim_var_cont[y, beta]
    y_covar <- sim_var_cont[y, covar]

    var.epsilon <- var.epsilon + (2*cov(y_beta*covars[,get(y_covar)], covars[,wes.betas])) + (2*cov(y_beta*covars[,get(y_covar)], covars[,array.betas])) + (2*cov(y_beta*covars[,get(y_covar)], covars[,cat_covar_2.betas]))

  }

  # r2 is the amount of variance explained by all simulated covariates for this phenotype
  r2 <- 0.05
  # So that makes the random factor (epsilon) of the model r2 - var.epsilon (scaled to between 0-1)
  var.epsilon <- var.epsilon * ((1-r2) / r2)

  covars[,sim_pheno:=rnorm(nrow(covars), mean = 0, sd = sqrt(var.epsilon))]
  covars[,sim_pheno:= sim_pheno +
    (PC1 * sim_var_cont[covar=='PC1', beta]) +
    (PC2 * sim_var_cont[covar=='PC2', beta]) +
    (PC3 * sim_var_cont[covar=='PC3', beta]) +
    (PC4 * sim_var_cont[covar=='PC4', beta]) +
    (PC5 * sim_var_cont[covar=='PC5', beta]) +
    (PC6 * sim_var_cont[covar=='PC6', beta]) +
    (PC7 * sim_var_cont[covar=='PC7', beta]) +
    (PC8 * sim_var_cont[covar=='PC8', beta]) +
    (PC9 * sim_var_cont[covar=='PC9', beta]) +
    (PC10 * sim_var_cont[covar=='PC10', beta]) +
    (age * sim_var_cont[covar=='age', beta]) +
    (age_squared * sim_var_cont[covar=='age_squared', beta]) +
    (sex * sim_var_cont[covar=='sex', beta]) +
    (quant_covar_1 * sim_var_cont[covar=='quant_covar_1', beta]) +
    wes.betas + array.betas + cat_covar_2.betas
    ]

  formatted_form <- as.formula(paste("sim_pheno", paste(c(paste0("PC",1:10), 'age', 'age_squared', 'sex', 'wes.batch', 'quant_covar_1', 'cat_covar_2', 'array'), collapse=" + "), sep = " ~ "))
  test_lm <- glance(glm(formatted_form, data = covars, family="gaussian"))
  cat(paste0('Variance explained by base covariates in simulated trait: ', 1 - (test_lm$deviance / test_lm$null.deviance)))

  # Now write three pieces of information:
  setkey(covars, 'eid')
  # 1. A base covariates file:
  cols <- c('eid', paste0("PC",1:40), 'age', 'sex', 'wes.batch', 'array', 'genetics_qc_pass')
  base_covar_file <- covars[,..cols]
  setnames(base_covar_file, cols, c('eid', paste('22009-0',1:40,sep='.'), '21003-0.0', '22001-0.0', 'wes.batch', '22000-0.0', 'genetics_qc_pass'))
  fwrite(base_covar_file, 'test_data/old_test_data/test_base_covariates.covariates', sep='\t', col.names = T, row.names = F, quote = F)

  # 2. Additional covariates file:
  cols <- c('FID','IID','quant_covar_1', 'cat_covar_2')
  fwrite(covars[,..cols], 'test_data/old_test_data/test_add_covariates.covariates', sep='\t', col.names = T, row.names = F, quote = F)

  # 3. A info file of covariate information including expected beta:
  setnames(wes_betas, c("wes.batch","wes.betas"), c("covar","beta"))
  setnames(array_betas, c("array","array.betas"), c("covar","beta"))
  setnames(cat_covar_betas, c("cat_covar_2","cat_covar_2.betas"), c("covar","beta"))
  covar_betas <- rbind(sim_var_cont[,c('covar','beta')], wes_betas, array_betas, cat_covar_betas)
  fwrite(covar_betas, 'test_data/old_test_data/covar_beta.info', sep='\t', col.names = T, row.names = F, quote = F)

  # And return the phenotype so we can continue to modify it with genetic data
  return(covars[,c('FID','IID','sim_pheno')])

}

simulate_genetic_effects <- function(phenotype) {

  phenotype[,sim_pheno_gt:=sim_pheno]
  genetic_variants <- fread('sim_data/sim_chromosome.bim')
  setnames(genetic_variants, names(genetic_variants), c('chrom','varID','cM','pos','REF','ALT'))
  effect_variants <- sample_n(genetic_variants, 1000)
  fwrite(effect_variants[,"varID"], 'sim_data/sim_chromosome.effectvars', row.names=F, col.names=F, quote=F, sep='\t')
  cmd <- "plink2 --bfile sim_data/sim_chromosome --extract sim_data/sim_chromosome.effectvars --out sim_data/sim_chromosome --export Av"
  system(cmd)

  effect_gts <- fread('sim_data/sim_chromosome.traw')
  setnames(effect_gts, paste(1000000:1009999, 1000000:1009999, sep='_'), as.character(1000000:1009999))

  # simulate random genetic effects for these variants
  effect_variants[,beta:=rnorm(nrow(effect_variants), mean=0, sd=0.15)]

  id_cols <- as.character(1000000:1009999)
  mod_pheno <- function(varID, beta) {
    curr_gts <- data.table::transpose(effect_gts[SNP == varID, ..id_cols], keep.names = 'IID')
    setnames(curr_gts, 'V1', 'gt')
    curr_gts[,gt:=2-gt]
    curr_gts[,effect:=gt * beta]
    phenotype[,sim_pheno_gt:=sim_pheno_gt + curr_gts[,effect]]
    test <- data.table(tidy(lm(sim_pheno_gt ~ gt, data=merge(phenotype, curr_gts, by ='IID'))))
    return(test[term=='gt',p.value])
  }
  effect_variants[,p:=mod_pheno(varID, beta),by=seq_len(nrow(effect_variants))]
  return(list(phenotype, effect_variants))
}

simulate_wes_data <- function(phenotype) {
  # Need to simulate the following in WES BGEN files for chromosomes 1 : 22 & X, with:
  #  * SNPs with specific annotations that vary based on
  #    + MAF
  #    + Consequence
  #    + Some reasonable number of additional variant categories
  #  * SNPs within the same GENE that add up to some effect
  #  * Independent SNPs (in the same GENE?) that add up to some effect
  #  * GENEs that add up to some effect

  phenotype[,sim_pheno_1_var:=sim_pheno_gt]
  phenotype[,sim_pheno_2_var:=sim_pheno_gt]
  bases <- c('A','T','C','G')
  aas <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","O","S","U","T","W","Y","V","B","Z","J")

  # Lets pick 40 causal genes randomly from the genome
  causal_genes <- copy(sample_n(transcripts, 16))
  # And then seperate them into Missense and LoF genes:
  causal_genes[,pheno_name:=c(rep('sim_pheno_1', 8),rep('sim_pheno_2',8))]
  causal_genes[,causality:=c(rep('MISS',4),rep('PTV',4), rep('MISS',4),rep('PTV',4))]
  causal_genes <- causal_genes[,c('ENST','causality','pheno_name')]

  gene_information_table <- data.table()

  # Need to write a header for each variant with these exactly columns:
  annote_header <- c("CHROM","POS","REF","ALT","varID","ogVarID","FILTER","AF","F_MISSING","AN","AC","MANE",
                     "ENST","ENSG","BIOTYPE","SYMBOL","CSQ","gnomAD_AF","CADD","REVEL","SIFT","POLYPHEN","LOFTEE",
                     "AA","AApos","PARSED_CSQ","MULTI","INDEL","MINOR","MAJOR","MAF","MAC")

  # Now we need to simulate variants in every gene we are going to look at.
  # I don't think it matters to the tools if they are actually in exons, just that they have the right notation
  # The number of variants that we can expect with MAF < 1e-3 follows a neg. bin. distribution with the following parameters:
  # (these were empirically determined from 470k UKBB data)
  sim_params <- data.table('CSQ'=c('MISS', 'PTV', 'SYN'), 'theta'=rep(1.172415, 3), 'mu'=c(72.31713,35.31489,163.3054))

  # I am going to simulate variants at ~10x the MAF they would actually be in the real data to ensure that the number
  # of variants in large enough to reasonably test the various pipelines. I don't think it should cause any issues...
  for (row in seq_len(nrow(transcripts))) {

    curr_gene <- transcripts[row]

    annote_table <- data.table('PARSED_CSQ' = c(rep("MISS", rnbinom(1, size=sim_params[CSQ == 'MISS', theta], mu=sim_params[CSQ == 'MISS', mu])),
                                                rep("PTV", rnbinom(1, size=sim_params[CSQ == 'PTV', theta], mu=sim_params[CSQ == 'PTV', mu])),
                                                rep("SYN", rnbinom(1, size=sim_params[CSQ == 'SYN', theta], mu=sim_params[CSQ == 'SYN', mu]))))

    # For non-causal genes, we downsample to make sure we have a small-ish test dataset that won't take a long time to run
    if (!curr_gene[,ENST] %in% causal_genes[,ENST]) {
      annote_table <- copy(sample_n(annote_table, nrow(annote_table) * 0.1))
      causal <- F
      causality <- NA
      phen <- NA
    } else {
      causal <- T
      causality <- causal_genes[ENST == curr_gene[,ENST],causality]
      phen <- causal_genes[ENST == curr_gene[,ENST],pheno_name]
    }

    if (nrow(annote_table) != 0) {

      annote_table[,CHROM:=paste0('chr',curr_gene[,chrom])]
      annote_table[,POS:=sample(curr_gene[,start]:curr_gene[,end], size = nrow(annote_table), replace = F)]
      annote_table[,REF:=sample(bases, size=nrow(annote_table), replace = T)]
      annote_table[,ALT:=sample(bases[!bases %in% REF], size=1, replace = T), by=seq_len(nrow(annote_table))]
      annote_table[,varID:=paste(CHROM, POS, REF, ALT, sep=':'), by=seq_len(nrow(annote_table))]
      annote_table[,ogVarID:=paste(CHROM, POS, REF, ALT, sep='_'), by=seq_len(nrow(annote_table))]
      annote_table[,FILTER:='PASS']
      annote_table[,AC:=sample(1:20, size = nrow(annote_table), replace = T, prob = 0.4*c(1:20)^-1.679)]  # -1.679 is from the empirical AC distribution
      annote_table[,AF:=AC / 20000]  # Always simulating 10k individuals
      annote_table[,F_MISSING:=0.123]  # Dummy variable to test functionality of logging in collapse variants, is not actually used
      annote_table[,AN:=20000]
      annote_table[,MANE:=curr_gene[,MANE]]
      annote_table[,ENST:=curr_gene[,ENST]]
      annote_table[,ENSG:=curr_gene[,ENSG]]
      annote_table[,BIOTYPE:=curr_gene[,BIOTYPE]]
      annote_table[,SYMBOL:=curr_gene[,SYMBOL]]
      annote_table[,CSQ:=paste(PARSED_CSQ, paste(paste0('testcsq', 1:sample(1:2, 1, replace=T)),collapse='&'), sep='&'), by=seq_len(nrow(annote_table))]
      annote_table[,gnomAD_AF:=0]
      annote_table[,CADD:=0.456]  # Dummy variable to test functionality of logging in collapse variants, is not actually used
      annote_table[,REVEL:=if_else(PARSED_CSQ=='MISS', 0.789, as.double(NA))]  # Dummy variable to test functionality of logging in collapse variants, is not actually used
      annote_table[,SIFT:=if_else(PARSED_CSQ=='MISS', 'tolerated_low_confidence(0.1)', as.character(NA))]  # Dummy variable to ensure pandas column type matching
      annote_table[,POLYPHEN:=if_else(PARSED_CSQ=='MISS', 'benign(0.05)', as.character(NA))]  # Dummy variable to ensure pandas column type matching
      annote_table[,LOFTEE:=if_else(PARSED_CSQ=='PTV','HC', as.character(NA))]
      annote_table[,AA:=if_else(PARSED_CSQ == 'MISS', paste(sample(aas, 2), collapse='/'), as.character(NA)), by=seq_len(nrow(annote_table))]
      annote_table[,AApos:=ifelse(PARSED_CSQ == 'MISS', POS - curr_gene[,start], as.integer(NA))]
      annote_table[,MULTI:=F]
      annote_table[,INDEL:=F]
      annote_table[,MINOR:=ALT]
      annote_table[,MAJOR:=REF]
      annote_table[,MAF:=AF]
      annote_table[,MAC:=AC]
      setkey(annote_table, 'POS')

      setcolorder(annote_table, annote_header)
      out_vep <- paste0('sim_data/chr',curr_gene[,chrom],'.filtered.vep.tsv')
      fwrite(annote_table, file=out_vep, append = file.exists(out_vep), quote=F, sep='\t', na='NA', row.names = F)

      out_traw <- paste0('sim_data/chr',curr_gene[,chrom],'.filtered.traw')

      # Ids formats weird...
      if (nrow(annote_table) == 1) {
        ids <- sample(phenotype[,IID], size=annote_table[,AC], replace=F)
      } else {
        ids <- unlist(sapply(annote_table[,AC], sample, x=phenotype[,IID], replace=F, simplify = F))
      }
      carriers <- data.table('CHR'=curr_gene[,chrom],
                             'SNP'=annote_table[,rep(varID,AC)],
                             'CM'=0,
                             'POS'=annote_table[,rep(POS,AC)],
                             'COUNTED'=annote_table[,rep(REF,AC)],
                             'ALT'=annote_table[,rep(ALT,AC)],
                             'IID'=ids,
                             'csq'=annote_table[,rep(PARSED_CSQ,AC)],
                             'gt'=1)
      csq_tots <- carriers[,sum(gt),by=c('csq','IID')]
      carriers_formatted <- data.table(pivot_wider(carriers, names_from = IID, values_from = gt, values_fill = 2, id_cols = c('CHR','SNP','CM','POS','COUNTED','ALT')))
      suppressWarnings(carriers_formatted <- carriers_formatted[,phenotype[!IID %in% carriers[,IID],IID]:=2])
      setcolorder(carriers_formatted,c('CHR','SNP','CM','POS','COUNTED','ALT',phenotype[,IID]))
      fwrite(carriers_formatted, file=out_traw, append = file.exists(out_traw), quote=F, sep='\t', na='NA', col.names=F, row.names = F)

      # Assign carriers and effect size (always + 2 / allele), if necessary
      if (curr_gene[,ENST] %in% causal_genes[,ENST]) {
        csq_count <- carriers[csq == causal_genes[ENST==curr_gene[,ENST],causality],sum(gt),by=IID]
        csq_count <- merge(phenotype[,'IID'], csq_count,all.x=T)
        curr_pheno <- paste0(causal_genes[ENST==curr_gene[,ENST],pheno_name],'_var')
        phenotype[,eval(curr_pheno):=get(curr_pheno)+csq_count[,if_else(is.na(V1),0, V1 * 2)]]
      }
    }
    current_gene_information <- data.table('ENST'=curr_gene[,ENST],
                                           'MISS'=nrow(annote_table[PARSED_CSQ == 'MISS']),
                                           'PTV'=nrow(annote_table[PARSED_CSQ == 'PTV']),
                                           'SYN'=nrow(annote_table[PARSED_CSQ == 'SYN']),
                                           'causal'=causal,'causality'=causality, 'pheno_name'=phen)
    if (nrow(annote_table) != 0) {
      current_gene_information[,'MISS_het':=nrow(csq_tots[csq=='MISS'])]
      current_gene_information[,'PTV_het':=nrow(csq_tots[csq=='PTV'])]
      current_gene_information[,'SYN_het':=nrow(csq_tots[csq=='SYN'])]
    } else {
      current_gene_information[,'MISS_het':=0]
      current_gene_information[,'PTV_het':=0]
      current_gene_information[,'SYN_het':=0]
    }
    gene_information_table <- rbind(gene_information_table, current_gene_information)
    if (row %% 100 == 0) {
      cat(paste0('Gene No ', row, '\n'))
    }

  }

  return(list(phenotype, gene_information_table))

}

# Make sure we start with actual genes:
transcripts <- fread('test_data/old_test_data/transcripts.tsv.gz')
transcripts <- transcripts[fail == F & chrom != 'Y']

# Simulate two independent phenotypes:
phenotype <- simulate_base_covariates_and_phenotype()
sim_gt <- simulate_genetic_effects(phenotype)
phenotype <- sim_gt[[1]]
effect_variants <- sim_gt[[2]]

sim_wes <- simulate_wes_data(phenotype)
phenotype <- sim_wes[[1]]
gene_info <- sim_wes[[2]]
fam <- data.table(FID=0, IID=1000000:1009999, col3=0, col4=0, col5=0, col6 =-9)

# Write all outputs
fwrite(phenotype[,c("FID","IID","sim_pheno_1_var","sim_pheno_2_var")], 'test_data/old_test_data/sim_pheno.txt', sep=' ', col.names = T, row.names = F, quote = F)
fwrite(gene_info, 'test_data/old_test_data/gene_beta.info', sep='\t', col.names = T, row.names = F, quote = F, na = 'NA')
fwrite(effect_variants, 'test_data/old_test_data/gt_beta.info', sep='\t', col.names = T, row.names = F, quote = F, na = 'NA')
fwrite(fam, 'sim_data/filtered.fam', sep=' ', col.names = F, row.names = F, quote = F)

# Read annotations back in so we can generate test data automatically
annotations <- data.table()
for (file in Sys.glob('sim_data/*.vep.tsv')) {
  annotations <- rbind(annotations, fread(file))
}

gene_info <- fread('test_data/old_test_data/gene_beta.info')

# Generate a json of testing parameters for filtering expressions:
expression_test_data <- data.table('variant_type' = c('PTV', 'MISS', 'SYN'), 'test_type' = 'expression', 'snp_list' = NA, 'gene_list' = NA)
get_test_data <- function(variant_type) {

  sub_table <- annotations[MAF < 0.001 & PARSED_CSQ == eval(variant_type)]
  expected_vars <- nrow(sub_table)
  chr_table <- sub_table[CHROM == 'chr1']
  test_gene <- data.table(chr_table[,table(ENST)])[N == max(N), ENST][1] # The [1] just ensures we only get one gene/variant in case of multiple fitting the criteria
  test_var <- chr_table[ENST == eval(test_gene)][MAC == max(MAC),varID][1]

  col <- paste(variant_type,'het',sep='_')
  gene_het_count <- gene_info[ENST == eval(test_gene),get(col)]
  var_het_count <- chr_table[ENST == eval(test_gene)][MAC == max(MAC),MAC]

  tot_gene_count <- length(chr_table[,table(ENST)])
  test_var_count <- nrow(chr_table[ENST == eval(test_gene)])
  tot_var_count <- nrow(chr_table)

  return(list(expected_vars, test_gene, test_var, gene_het_count, var_het_count, tot_gene_count, test_var_count, tot_var_count))

}
expression_test_data[,c("expected_vars", "test_gene", "test_var", "gene_het_count", "var_het_count", "tot_gene_count", "test_var_count","tot_var_count"):=get_test_data(variant_type),by=seq_len(nrow(expression_test_data))]

# Generate a json of testing parameters for SNP and GENE lists:
list_test_data <- data.table('variant_type' = c('PTV', 'MISS', 'PTV', 'MISS'),
                             'test_type' = c('gene', 'gene', 'snp', 'snp'),
                             'snp_list' = c(NA, NA, 'snp_list.pheno_1_PTV.txt', 'snp_list.pheno_1_MISS.txt'),
                             'gene_list' = c('gene_list.pheno_1_PTV.txt', 'gene_list.pheno_1_MISS.txt', NA, NA))

# This function calculates numbers of variants when collapsing on gene / snp masks
get_het_count <- function(genes, variant_type, maf_filter) {

  snp_table <- data.table()
  for (gene in genes) {
    curr_traw <- fread(paste0('sim_data/chr',transcripts[ENST == gene,chrom],'.filtered.traw'),col.names = c('CHROM','SNP','cM','POS','REF','ALT',1000000:1009999))
    curr_snps <- annotations[ENST == gene & PARSED_CSQ == variant_type & MAF < maf_filter,varID]
    cols <- c('SNP',1000000:1009999)
    curr_traw <- data.table(pivot_longer(curr_traw[SNP %in% curr_snps, ..cols], -SNP, values_to="GT", names_to="IID"))
    curr_traw <- curr_traw[GT != 2]
    snp_table <- rbind(snp_table, curr_traw)
  }
  return(nrow(snp_table[,sum(GT),by="IID"]))
}

get_test_data <- function(variant_type, snp_list, gene_list) {

  if (!is.na(snp_list)) {
    sub_table <- annotations[PARSED_CSQ == eval(variant_type)]
    filter_list <- annotations[ENST %in% gene_info[causality == eval(variant_type) & pheno_name == 'sim_pheno_1',ENST] & PARSED_CSQ == eval(variant_type),varID]
    sub_table <- sub_table[varID %in% filter_list]
    test_gene <- 'ENST00000000000'
    gene_het_count <- get_het_count(sub_table[,unique(ENST)], variant_type, 1)
  } else if (!is.na(gene_list)) {
    sub_table <- annotations[MAF < 0.001 & PARSED_CSQ == eval(variant_type)]
    filter_list <- gene_info[causal == T & causality == eval(variant_type) & pheno_name == 'sim_pheno_1',ENST]
    sub_table <- sub_table[ENST %in% filter_list]
    test_gene <- 'ENST99999999999'
    gene_het_count <- get_het_count(sub_table[,unique(ENST)], variant_type, 0.001)
  }

  expected_vars <- nrow(sub_table)
  test_var <- sub_table[MAC == max(MAC),varID][1]
  var_het_count <- sub_table[varID == eval(test_var),MAC]
  tot_gene_count <- 1  # Always 1 because we collapse down to a single gene (ENST00000000000 / ENST99999999999)
  test_var_count <- nrow(sub_table)
  tot_var_count <- nrow(sub_table)

  return(list(expected_vars, test_gene, test_var, gene_het_count, var_het_count, tot_gene_count, test_var_count, tot_var_count))

}
list_test_data[,c("expected_vars", "test_gene", "test_var", "gene_het_count", "var_het_count", "tot_gene_count", "test_var_count","tot_var_count"):=get_test_data(variant_type, snp_list, gene_list),by=seq_len(nrow(list_test_data))]

# And cat together and print
cat(toJSON(rbind(expression_test_data, list_test_data), pretty=T, na='null'), file= "test_data/old_test_data/expression_test_data.json")

# Generate a set of gene / SNP lists for testing gene / SNP collapsing:
fwrite(annotations[ENST %in% gene_info[causal == T & causality == 'MISS' & pheno_name == 'sim_pheno_1',ENST] & PARSED_CSQ == 'MISS','varID'], "test_data/old_test_data/snp_list.pheno_1_MISS.txt", row.names=F, col.names=F, quote=F, sep="\t")
fwrite(annotations[ENST %in% gene_info[causal == T & causality == 'PTV' & pheno_name == 'sim_pheno_1',ENST] & PARSED_CSQ == 'PTV','varID'], "test_data/old_test_data/snp_list.pheno_1_PTV.txt", row.names=F, col.names=F, quote=F, sep="\t")
fwrite(annotations[ENST %in% gene_info[causal == T & causality == 'MISS' & pheno_name == 'sim_pheno_2',ENST] & PARSED_CSQ == 'MISS','varID'], "test_data/old_test_data/snp_list.pheno_2_MISS.txt", row.names=F, col.names=F, quote=F, sep="\t")
fwrite(annotations[ENST %in% gene_info[causal == T & causality == 'PTV' & pheno_name == 'sim_pheno_2',ENST] & PARSED_CSQ == 'PTV','varID'], "test_data/old_test_data/snp_list.pheno_2_PTV.txt", row.names=F, col.names=F, quote=F, sep="\t")

fwrite(transcripts[ENST %in% gene_info[causal == T & causality == 'MISS' & pheno_name == 'sim_pheno_1',ENST],'SYMBOL'], "test_data/old_test_data/gene_list.pheno_1_MISS.txt", row.names=F, col.names=F, quote=F, sep="\t")
fwrite(transcripts[ENST %in% gene_info[causal == T & causality == 'PTV' & pheno_name == 'sim_pheno_1',ENST],'SYMBOL'], "test_data/old_test_data/gene_list.pheno_1_PTV.txt", row.names=F, col.names=F, quote=F, sep="\t")
fwrite(transcripts[ENST %in% gene_info[causal == T & causality == 'MISS' & pheno_name == 'sim_pheno_2',ENST],'SYMBOL'], "test_data/old_test_data/gene_list.pheno_2_MISS.txt", row.names=F, col.names=F, quote=F, sep="\t")
fwrite(transcripts[ENST %in% gene_info[causal == T & causality == 'PTV' & pheno_name == 'sim_pheno_2',ENST],'SYMBOL'], "test_data/old_test_data/gene_list.pheno_2_PTV.txt", row.names=F, col.names=F, quote=F, sep="\t")
