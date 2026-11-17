# Figure 2: Root and shoot nitrogen content and C:N ratio at day 14
# Four panels showing mean ± SE for Basalt, Oakley and F1 under HN and LN.
# Tukey letters compare genotypes within each treatment.
# Brackets show pairwise genotype comparisons under LN only.
# Input: BIGS_C_N_balanced_3rep(1).xlsx

options(stringsAsFactors = FALSE)

# Load packages
for (pkg in c('readxl','multcompView')) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos='https://cloud.r-project.org')
  }
  if (!requireNamespace(pkg, quietly = TRUE)) stop('Cannot load ', pkg)
}

# Read and check the workbook
input_file <- file.choose()
dat <- as.data.frame(readxl::read_excel(input_file, sheet=1))
required <- c('Genotype','N_treatment','Organ','N_pct','C_pct','C_to_N')
if (!all(required %in% names(dat))) stop('Expected columns: ', paste(required,collapse=', '),
                                         '\nFound: ', paste(names(dat),collapse=', '))
for (v in required[4:6]) dat[[v]] <- as.numeric(dat[[v]])
dat$Genotype <- trimws(as.character(dat$Genotype))
dat$N_treatment <- trimws(as.character(dat$N_treatment))
dat$Organ <- trimws(as.character(dat$Organ))
dat <- dat[dat$Genotype %in% c('Basalt','Oakley','Basalt x Oakley') &
             dat$N_treatment %in% c('High N','Low N') &
             dat$Organ %in% c('Root','Shoot'), ,drop=FALSE]
if (nrow(dat) == 0 || anyNA(dat[,required])) stop('Missing/invalid data; check the workbook.')
genotypes <- c('Basalt','Oakley','Basalt x Oakley')
nitrogens <- c('High N','Low N')
organs <- c('Root','Shoot')
for (g in genotypes) for (n in nitrogens) for (o in organs) {
  count <- sum(dat$Genotype==g & dat$N_treatment==n & dat$Organ==o)
  if (count != 3) stop('Expected 3 replicates for ',g,', ',n,', ',o,'; found ',count)
}

out_dir <- file.path(dirname(input_file),'Figure2_CN_FINAL_4PANELS')
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Compact letters from Tukey-adjusted pairwise p-values
get_letters <- function(pvals) {
  letter_values <- multcompView::multcompLetters(pvals, threshold=0.05, Letters=base::letters)$Letters
  setNames(as.character(letter_values), names(letter_values))
}

specs <- list(
  list(organ='Root',  variable='N_pct',  ylabel='Root nitrogen (%)'),
  list(organ='Root',  variable='C_to_N', ylabel='Root C:N ratio'),
  list(organ='Shoot', variable='N_pct',  ylabel='Shoot nitrogen (%)'),
  list(organ='Shoot', variable='C_to_N', ylabel='Shoot C:N ratio')
)

# Calculate means, standard errors, ANOVA and Tukey comparisons
summary_rows <- list(); comparison_rows <- list(); anova_rows <- list()
results <- vector('list',4)
for (k in seq_along(specs)) {
  sp <- specs[[k]]
  sub <- dat[dat$Organ==sp$organ, ,drop=FALSE]
  sub$Value <- sub[[sp$variable]]
  sub$Genotype <- factor(sub$Genotype,levels=genotypes)
  sub$N_treatment <- factor(sub$N_treatment,levels=nitrogens)
  # The two-factor ANOVA is saved separately; plot annotations use within-treatment tests.
  fit2 <- aov(Value ~ Genotype * N_treatment, data=sub)
  at <- as.data.frame(summary(fit2)[[1]])
  at$Effect <- rownames(at); at$Panel <- sp$ylabel
  anova_rows[[k]] <- at
  means <- matrix(NA_real_,3,2,dimnames=list(genotypes,nitrogens))
  ses <- means
  labs <- matrix('',3,2,dimnames=list(genotypes,nitrogens))
  pLN <- setNames(rep(NA_real_,3),c('Basalt-Oakley','Basalt-F1','Oakley-F1'))
  for (j in seq_along(nitrogens)) {
    nt <- nitrogens[j]
    d <- sub[sub$N_treatment==nt,,drop=FALSE]
    d$Genotype <- droplevels(d$Genotype)
    for (i in seq_along(genotypes)) {
      v <- d$Value[d$Genotype==genotypes[i]]
      means[i,j] <- mean(v)
      ses[i,j] <- sd(v)/sqrt(length(v))
      summary_rows[[length(summary_rows)+1]] <- data.frame(Panel=sp$ylabel,Organ=sp$organ,
        Trait=sp$variable,Genotype=genotypes[i],Nitrogen=nt,n=length(v),Mean=means[i,j],SE=ses[i,j])
    }
    fit <- aov(Value ~ Genotype,data=d)
    tuk <- TukeyHSD(fit,'Genotype')$Genotype
    ps <- tuk[,'p adj']; names(ps) <- rownames(tuk)
    lettermap <- get_letters(ps)
    labs[,j] <- unname(lettermap[genotypes])
    if (anyNA(labs[,j])) stop('Letter assignment failed for ',sp$ylabel,' ',nt)
    for (ii in seq_len(nrow(tuk))) {
      comparison_rows[[length(comparison_rows)+1]] <- data.frame(Panel=sp$ylabel,
        Nitrogen=nt,Comparison=rownames(tuk)[ii],Difference=tuk[ii,'diff'],
        LowerCI=tuk[ii,'lwr'],UpperCI=tuk[ii,'upr'],Tukey_p=tuk[ii,'p adj'])
    }
    if (nt=='Low N') {
      # TukeyHSD reports each difference as the second level minus the first.
      pLN['Basalt-Oakley'] <- unname(ps['Oakley-Basalt'])
      pLN['Basalt-F1'] <- unname(ps['Basalt x Oakley-Basalt'])
      pLN['Oakley-F1'] <- unname(ps['Basalt x Oakley-Oakley'])
      if (anyNA(pLN)) stop('Missing LN pairwise comparisons for ',sp$ylabel)
    }
  }
  results[[k]] <- list(spec=sp,means=means,ses=ses,letters=labs,pLN=pLN)
}

# Save statistical tables
write.csv(do.call(rbind,summary_rows),file.path(out_dir,'Descriptive_means_SE.csv'),row.names=FALSE)
write.csv(do.call(rbind,comparison_rows),file.path(out_dir,'Within_N_Tukey_comparisons.csv'),row.names=FALSE)
write.csv(do.call(rbind,anova_rows),file.path(out_dir,'Two_way_ANOVA.csv'),row.names=FALSE)

# Plot colours and significance labels
HN <- '#E58D7D'; LN <- '#59B5AA'
star <- function(p) if (is.na(p)) 'NA' else if (p < 0.001) '***' else
  if (p < 0.01) '**' else if (p < 0.05) '*' else 'ns'

plot_panel <- function(res,tag) {
  m <- res$means; se <- res$ses; lab <- res$letters
  ymax <- max(m+se)
  upper <- ymax*1.38
  par(mar=c(4.7,5.5,2.9,1.3), xpd=FALSE)
  plot(NA,xlim=c(0.45,3.55),ylim=c(0,upper),xaxt='n',xlab='',ylab='',
       bty='l',las=1,cex.axis=0.88)
  axis(1,at=1:3,labels=c('Basalt','Oakley',expression(F[1])),tick=FALSE,
       line=0.25,cex.axis=0.95)
  title(ylab=res$spec$ylabel,cex.lab=1.02,line=3.2)
  width <- 0.31; gap <- 0.02
  xpos <- cbind(1:3-(width+gap)/2,1:3+(width+gap)/2)
  for (i in 1:3) for (j in 1:2) {
    x <- xpos[i,j]; y <- m[i,j]; err <- se[i,j]
    rect(x-width/2,0,x+width/2,y,col=c(HN,LN)[j],border='#555555',lwd=0.85)
    arrows(x,y-err,x,y+err,angle=90,code=3,length=0.045,lwd=0.85)
    text(x,y+err+0.025*upper,labels=lab[i,j],cex=0.87,font=2)
  }
  # Three LN genotype comparisons are shown above the bars.
  br0 <- max(m+se)+0.105*upper
  increment <- 0.077*upper
  brackets <- list(c(1,2,'Basalt-Oakley'),c(2,3,'Oakley-F1'),c(1,3,'Basalt-F1'))
  for (q in seq_along(brackets)) {
    b <- brackets[[q]]; x1 <- as.numeric(b[1]); x2 <- as.numeric(b[2]);
    yy <- br0+(q-1)*increment
    segments(x1,yy-0.016*upper,x1,yy,lwd=0.8)
    segments(x1,yy,x2,yy,lwd=0.8)
    segments(x2,yy,x2,yy-0.016*upper,lwd=0.8)
    text((x1+x2)/2,yy+0.013*upper,star(res$pLN[b[3]]),cex=0.8)
  }
  mtext(paste0('(', tag, ') ', res$spec$ylabel),
        side=3, adj=0, line=1.25, font=2, cex=0.95)
}

# Arrange the four panels with a shared legend
render <- function() {
  old <- par(no.readonly=TRUE)
  on.exit(par(old), add=TRUE)
  layout(matrix(c(1,1,2,3,4,5), nrow=3, byrow=TRUE),
         heights=c(0.18,1,1), widths=c(1,1))
  par(family='sans', cex=1.10, oma=c(0,0,0,0))
  par(mar=c(0,0,0,0), xpd=NA)
  plot.new()
  legend('center', legend=c('High N','Low N'), fill=c(HN,LN),
         border='#555555', horiz=TRUE, bty='n', cex=1.0,
         x.intersp=0.8, y.intersp=1)
  for (k in seq_len(4)) plot_panel(results[[k]], base::letters[k])
}

# Export the figure
png_file <- file.path(out_dir,'Figure2_CN_Four_Panels_600dpi.png')
pdf_file <- file.path(out_dir,'Figure2_CN_Four_Panels_Vector.pdf')
png(png_file,width=13,height=10,units='in',res=600,type='cairo',bg='white')
tryCatch(render(),finally=dev.off())
pdf(pdf_file,width=13,height=10,onefile=FALSE,useDingbats=FALSE)
tryCatch(render(),finally=dev.off())
cat('\nFiles saved in:\n',out_dir,'\n',png_file,'\n',pdf_file,'\n')
