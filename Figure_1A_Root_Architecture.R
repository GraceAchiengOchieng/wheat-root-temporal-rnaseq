# Figure 1A: Root architecture at day 14
# Root length, surface area, volume and diameter under HN and LN
# Bars show mean ± SE. Letters indicate Tukey groups within each treatment.
# Brackets show pairwise genotype comparisons under LN, including ns.
# Inputs: Group_means_SD_SE.csv, Tukey_within_treatment.csv,
#         Tukey_letters_within_treatment.csv

# Packages
pkgs <- c('ggplot2','dplyr','readr','patchwork')
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

# Select the three input files
read_checked <- function(prompt, required) {
  repeat {
    message('\nSELECT: ', prompt)
    path <- file.choose()
    d <- readr::read_csv(path, show_col_types=FALSE, name_repair='minimal')
    names(d) <- trimws(sub('^\\ufeff','',names(d)))
    if (ncol(d)==1 && grepl(';', names(d)[1], fixed=TRUE)) {
      d <- readr::read_delim(path, delim=';', show_col_types=FALSE)
      names(d) <- trimws(names(d))
    }
    if (all(required %in% names(d))) return(d)
    message('Wrong CSV. Missing: ', paste(setdiff(required,names(d)),collapse=', '), '\nTry again.')
  }
}
means <- read_checked('Group_means_SD_SE.csv', c('Trait','Genotype','Treatment','mean','se'))
tukey <- read_checked('Tukey_within_treatment.csv', c('Trait','Treatment','contrast','p.value'))
letters <- read_checked('Tukey_letters_within_treatment.csv', c('Trait','Genotype','Treatment','Letter'))

# Traits and plot labels
traits <- c('Root length (cm)','Surface area (cm²)','Root volume (cm³)','Root diameter (mm)')
titles <- c('Total root length','Root surface area','Root volume','Root diameter')
ylabels <- c('Root length (cm)','Surface area (cm²)','Root volume (cm³)','Root diameter (mm)')
for (name in c('means','tukey','letters')) {
  obj <- get(name)
  obj$Trait <- trimws(as.character(obj$Trait))
  if (!all(traits %in% obj$Trait)) stop('Missing expected traits in ',name, ': ',paste(setdiff(traits,obj$Trait),collapse='; '))
  assign(name,obj)
}
means <- means |> dplyr::filter(Trait %in% traits) |>
  dplyr::mutate(Genotype=as.character(Genotype),Treatment=as.character(Treatment))
letters <- letters |> dplyr::filter(Trait %in% traits) |>
  dplyr::mutate(Genotype=as.character(Genotype),Treatment=as.character(Treatment))
stopifnot(nrow(means)==24, nrow(letters)==24)
if (anyDuplicated(means[c('Trait','Genotype','Treatment')]) ||
    anyDuplicated(letters[c('Trait','Genotype','Treatment')])) stop('Duplicated group keys in means or letters.')
means <- dplyr::left_join(means,letters,by=c('Trait','Genotype','Treatment'))
if (anyNA(means$Letter)) stop('Missing Tukey letters for some groups.')
means$Genotype <- factor(means$Genotype,levels=c('Basalt','Oakley','F1'))
means$Treatment <- factor(means$Treatment,levels=c('HN','LN'))
if (anyNA(means$Genotype) || anyNA(means$Treatment)) stop('Unexpected genotype or treatment labels.')

# Tukey comparisons between genotypes under LN
wanted <- c('Basalt - Oakley','Basalt - F1','Oakley - F1')
ln <- tukey |> dplyr::filter(Trait %in% traits,Treatment=='LN',contrast %in% wanted)
if (nrow(ln)!=12 || anyDuplicated(ln[c('Trait','contrast')])) stop('Expected 3 unique LN Tukey contrasts per trait; inspect Tukey CSV.')
ln$label <- ifelse(ln$p.value<.001,'***',ifelse(ln$p.value<.01,'**',ifelse(ln$p.value<.05,'*','ns')))
# The input p-values are already Tukey-adjusted.

# Save the tables used for plotting
out <- file.path(getwd(),'Figure_1A_FINAL')
dir.create(out,showWarnings=FALSE,recursive=TRUE)
readr::write_csv(ln,file.path(out,'LN_genotype_Tukey_comparisons.csv'))
readr::write_csv(means,file.path(out,'Plotted_group_means_SE_and_letters.csv'))

pal <- c(HN='#D98E9A',LN='#4AA9A6')

# Plot each root trait
# Dodged bar centres are shifted by 0.195 from each genotype position.
make_panel <- function(i) {
  trait <- traits[i]
  d <- dplyr::filter(means,Trait==trait)
  t <- dplyr::filter(ln,Trait==trait)
  top <- max(d$mean+d$se)
  span <- max(top,0.1)
  d$letter_y <- d$mean+d$se+span*.035
  t$level <- match(t$contrast,c('Basalt - Oakley','Oakley - F1','Basalt - F1'))
  t$y <- top+span*(.14+.115*(t$level-1))
  t$x1 <- c('Basalt - Oakley'=1,'Basalt - F1'=1,'Oakley - F1'=2)[t$contrast]+.195
  t$x2 <- c('Basalt - Oakley'=2,'Basalt - F1'=3,'Oakley - F1'=3)[t$contrast]+.195
  tick <- span*.018
  ymax <- max(t$y)+span*.13

  ggplot2::ggplot(d,ggplot2::aes(x=Genotype,y=mean,fill=Treatment))+
    ggplot2::geom_col(position=ggplot2::position_dodge(width=.78),width=.72,
                      colour='#333333',linewidth=.3)+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-se,ymax=mean+se),
                           position=ggplot2::position_dodge(width=.78),width=.12,linewidth=.5)+
    ggplot2::geom_text(ggplot2::aes(y=letter_y,label=Letter,group=Treatment),
                       position=ggplot2::position_dodge(width=.78),size=4.6,
                       vjust=0,fontface='plain',show.legend=FALSE)+
    ggplot2::geom_segment(data=t,ggplot2::aes(x=x1,xend=x2,y=y,yend=y),inherit.aes=FALSE,linewidth=.48)+
    ggplot2::geom_segment(data=t,ggplot2::aes(x=x1,xend=x1,y=y-tick,yend=y),inherit.aes=FALSE,linewidth=.48)+
    ggplot2::geom_segment(data=t,ggplot2::aes(x=x2,xend=x2,y=y-tick,yend=y),inherit.aes=FALSE,linewidth=.48)+
    ggplot2::geom_text(data=t,ggplot2::aes(x=(x1+x2)/2,y=y+span*.014,label=label),
                       inherit.aes=FALSE,size=4.4,vjust=0)+
    ggplot2::scale_fill_manual(values=pal,breaks=c('HN','LN'))+
    ggplot2::scale_x_discrete(labels=c('Basalt','Oakley',expression(F[1])))+
    ggplot2::scale_y_continuous(limits=c(0,ymax),expand=ggplot2::expansion(mult=c(0,0)))+
    ggplot2::labs(title=paste0('(',panel_ids[i],') ',titles[i]),x=NULL,y=ylabels[i],fill='Treatment')+
    ggplot2::theme_classic(base_size=14)+
    ggplot2::theme(plot.title=ggplot2::element_text(size=15,face='bold',hjust=0,margin=ggplot2::margin(b=8)),
                   axis.title.y=ggplot2::element_text(size=14),
                   axis.text=ggplot2::element_text(size=13,colour='black'),
                   axis.line=ggplot2::element_line(linewidth=.55),
                   legend.position='bottom',legend.title=ggplot2::element_text(size=13),
                   legend.text=ggplot2::element_text(size=13),
                   plot.margin=ggplot2::margin(8,16,8,8))
}

# Combine the four panels
panel_ids <- c('a','b','c','d')
panels <- lapply(seq_along(traits),make_panel)
fig <- (panels[[1]]+panels[[2]])/(panels[[3]]+panels[[4]])+
  patchwork::plot_layout(guides='collect') &
  ggplot2::theme(legend.position='bottom')
print(fig)

# Export figure
ggplot2::ggsave(file.path(out,'Figure_1A_FINAL_600dpi.png'),fig,width=14,height=10.7,units='in',dpi=600,bg='white')
ggplot2::ggsave(file.path(out,'Figure_1A_FINAL_vector.pdf'),fig,width=14,height=10.7,units='in',device=grDevices::cairo_pdf,bg='white')

message('\nSaved to: ',normalizePath(out))
message('Letters show genotype comparisons within each treatment.')
message('Brackets show Tukey-adjusted genotype comparisons under LN, including ns.')
message('HN vs LN comparisons and individual replicate values are not shown.')
