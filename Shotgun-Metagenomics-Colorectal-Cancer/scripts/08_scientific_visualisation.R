# PROJECT 4 | Scientific figure atlas
# Run from the project root with: Rscript scripts/08_scientific_visualisation.R
# Figures are generated from local project data and saved analysis results.
# Outputs: figures/reproduced/*.webp and figures/reproduced/figure_manifest.csv.
# Source-project accession is used for provenance. Study-level independence is not presumed.

options(stringsAsFactors = FALSE)
ROOT <- normalizePath(Sys.getenv('PROJECT4_ROOT', unset = getwd()), winslash = '/', mustWork = TRUE)
OUT <- file.path(ROOT, 'figures', 'reproduced')
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(20261027)

required_packages <- c('data.table','ggplot2','scales','ragg','vegan','patchwork','ComplexHeatmap','circlize','RColorBrewer')
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop('Missing R packages: ', paste(missing_packages, collapse = ', '),
       '. Install dependencies before running this script.')
}
library(data.table); library(ggplot2); library(scales); library(vegan)
library(ComplexHeatmap); library(circlize); library(grid)
ht_opt$message <- FALSE

PALETTE <- c(healthy='#087F8C',adenoma='#E3A546',CRC='#AC3656')
POS <- '#B54458'; NEG <- '#137D91'; DARK <- '#181C27'
manifest <- data.table(Figure=character(),Title=character(),Status=character(),File=character(),Note=character())
record <- function(id,title,status,path='',note='') manifest <<- rbind(manifest,data.table(Figure=id,Title=title,Status=status,File=path,Note=note))
KEEP_FIGURES <- sprintf('Figure_%02d', c(1:6, 8:13, 15, 17:24))
requested_figures <- Sys.getenv('PROJECT4_FIGURES', unset='')
if(nzchar(requested_figures)) KEEP_FIGURES <- trimws(strsplit(requested_figures, ',', fixed=TRUE)[[1]])
with_fig <- function(id,title,fun) {
  if (!id %in% KEEP_FIGURES) return(invisible(NULL))
  tryCatch({p <- fun(); record(id,title,'CREATED',p,'')},error=function(e){record(id,title,'SKIPPED','',conditionMessage(e))})
}
base_theme <- function(sz=14) theme_minimal(base_size=sz)+theme(
  text=element_text(colour='black',face='bold'),
  axis.title=element_text(colour='black',face='bold'),axis.text=element_text(colour='black',face='bold'),
  plot.title=element_text(face='bold',colour='black',size=sz+3),
  plot.subtitle=element_blank(),plot.caption=element_blank(),
  legend.title=element_text(face='bold',colour='black'),legend.text=element_text(face='bold',colour='black'),
  panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),
  plot.background=element_rect(fill='white',colour=NA),panel.background=element_rect(fill='white',colour=NA),
  strip.text=element_text(face='bold',colour='black'),plot.title.position='plot')
save_plot <- function(id,plot,w=12,h=8) {
  path <- file.path(OUT,paste0(id,'.webp'))
  ggplot2::ggsave(path,plot,width=w,height=h,dpi=220,device=ragg::agg_webp,bg='white',limitsize=FALSE)
  path
}
heat_col <- circlize::colorRamp2(c(0,0.5,1),c('#F5F7FA','#9EC8DB','#163E6A'))
heat_save <- function(id,mat,anno,header,w=15,h=8,legend='Log10 relative abundance',breaks=NULL,annotation_legend=TRUE,legend_side='right') {
  stopifnot(ncol(mat)==nrow(anno),identical(colnames(mat),anno$sample))
  q <- if (is.null(breaks)) as.numeric(quantile(mat,c(0,.50,.98),na.rm=TRUE)) else breaks
  if (length(unique(q)) < 3L) q <- c(min(mat),mean(range(mat)),max(mat))
  if(length(unique(q))<3L) stop('Heatmap contains no numerical variation')
  cfun <- circlize::colorRamp2(q,c('#F4F7FA','#91BFD4','#153957'))
  countries <- as.character(anno$country)
  countries[is.na(countries) | !nzchar(countries)] <- 'Not recorded'
  country_levels <- sort(unique(countries))
  country_cols <- setNames(grDevices::hcl.colors(length(country_levels),'Dark 3'),country_levels)
  # Annotation colours are keyed to source project but not treated as biological replicates.
  studies <- sort(unique(as.character(anno$study_accession)))
  study_cols <- setNames(grDevices::hcl.colors(length(studies),'Set 2'),studies)
  top <- ComplexHeatmap::HeatmapAnnotation(
    Clinical_group=as.character(anno$group),
    Country=countries,
    Source_project=as.character(anno$study_accession),
    col=list(Clinical_group=PALETTE,Country=country_cols,Source_project=study_cols),
    show_legend=annotation_legend,
    show_annotation_name=TRUE,
    annotation_name_gp=grid::gpar(col='black',fontface='bold',fontsize=10),
    simple_anno_size=grid::unit(3.5,'mm'))
  hm <- ComplexHeatmap::Heatmap(mat,name=legend,col=cfun,
    top_annotation=top,show_column_names=FALSE,show_row_names=TRUE,
    row_names_gp=grid::gpar(col='black',fontface='bold',fontsize=10),
    row_names_max_width=grid::unit(9,'cm'),
    column_split=factor(as.character(anno$group),levels=c('healthy','adenoma','CRC')),
    cluster_columns=FALSE,cluster_rows=TRUE,show_row_dend=FALSE,
    column_title_gp=grid::gpar(col='black',fontface='bold',fontsize=12),
    row_title_gp=grid::gpar(col='black',fontface='bold'),
    heatmap_legend_param=list(title_gp=grid::gpar(col='black',fontface='bold'),
      labels_gp=grid::gpar(col='black',fontface='bold')),
    use_raster=ncol(mat)>120,raster_quality=3,border=FALSE)
  file <- file.path(OUT,paste0(id,'.webp'))
  ragg::agg_webp(file,width=w,height=h,units='in',res=220,background='white')
  on.exit(grDevices::dev.off(),add=TRUE)
  # ComplexHeatmap::draw manages its own top-level viewport.
  # A nested pushViewport/popViewport around draw causes the reported grid error.
  ComplexHeatmap::draw(hm,heatmap_legend_side=legend_side,annotation_legend_side='right',
    merge_legends=TRUE,padding=grid::unit(c(9,6,6,6),'mm'))
  invisible(file)
}

# Read the original eligible metadata and numeric relative-signal profiles.
mfile <- file.path(ROOT,'results','stage1','metadata_eligible_primary.csv')
sfile <- file.path(ROOT,'data','species_signal_2340_CRC_cohort_20240322.tab')
if (!file.exists(mfile)||!file.exists(sfile)) stop('Original discovery inputs not found in the configured ROOT: ',ROOT)
md0 <- fread(mfile); sp <- fread(sfile,showProgress=FALSE)
stopifnot(all(c('sample','class','study_accession')%in%names(md0)), 'msp_id'%in%names(sp))
md0 <- md0[class%in%c('healthy','adenoma','CRC') & !is.na(study_accession)]
stopifnot(!anyDuplicated(md0$sample),all(md0$sample%in%names(sp)),!anyDuplicated(sp$msp_id))
md0[,group:=factor(class,levels=c('healthy','adenoma','CRC'))]
md0[,study_accession:=as.character(study_accession)]
setorder(md0,study_accession,group,sample)
input_sample_ids <- md0$sample
raw <- t(as.matrix(sp[, ..input_sample_ids])); storage.mode(raw)<-'double'
rownames(raw)<-md0$sample;colnames(raw)<-sp$msp_id
stopifnot(all(is.finite(raw)),all(raw>=0),all(rowSums(raw)>0))
relative <- sweep(raw,1,rowSums(raw),'/')
rm(raw,sp)

both <- md0[,.(crc=sum(group=='CRC'),healthy=sum(group=='healthy')),by=study_accession][crc>0&healthy>0,study_accession]
core <- md0[study_accession%in%both & group%in%c('CRC','healthy')]
core_idx <- match(core$sample,rownames(relative))
stopifnot(!anyNA(core_idx))
A <- relative[core_idx,,drop=FALSE]
feature_prevalence <- colMeans(A>0)

# ------------------------------------------------------------
# Taxonomy mapping from saved project annotations.
# Searches small annotation/Stage 4 tables for an exact msp_id -> species pair.
# Species figures are skipped when no mapped species names are available.
# A preferred manual file is results/stage4/msp_taxonomy.csv containing
# msp_id,species_name. The numeric profiles alone do NOT encode species names.
# ------------------------------------------------------------
get_taxonomy <- function() {
  explicit <- file.path(ROOT,'results','stage4','msp_taxonomy.csv')
  places <- c(file.path(ROOT,'results','stage4'),file.path(ROOT,'data'))
  candidates <- unique(c(explicit,unlist(lapply(places,function(d) {
    if(!dir.exists(d))return(character())
    list.files(d,pattern='(taxon|taxonom|species|msp|annotat).*\\.(csv|tsv|tab|txt)$',
               ignore.case=TRUE,recursive=TRUE,full.names=TRUE)
  }))))
  candidates <- candidates[file.exists(candidates)]
  # Do not read the enormous 2,340-column abundance matrix as taxonomy.
  candidates <- candidates[normalizePath(candidates,winslash='/',mustWork=FALSE)!=
    normalizePath(sfile,winslash='/',mustWork=FALSE)]
  accepted <- list(); audit <- data.table(File=character(),Mapped=integer())
  for(f in candidates) {
    if(file.info(f)$size>40e6)next
    d <- tryCatch(data.table::fread(f,showProgress=FALSE),error=function(e)NULL)
    if(is.null(d) || !nrow(d))next
    nn <- names(d); low<-tolower(nn)
    id <- which(low %in% c('msp_id','msp','mspid','feature_id','species_id'))
    lab <- which(low %in% c('species_name','species','scientific_name','taxon_name','species_label'))
    if(!length(id) || !length(lab))next
    x<-data.table(msp_id=as.character(d[[id[1]]]),
                  species_name=trimws(as.character(d[[lab[1]]])))
    x<-x[msp_id %in% colnames(A) & !is.na(species_name) & nzchar(species_name) &
      !grepl('^(na|unknown|unclassified|uncultured|msp[_ ]?[0-9]+)',species_name,ignore.case=TRUE)]
    x<-unique(x,by=c('msp_id','species_name'))
    if(!nrow(x))next
    x[,Source:=f]
    accepted[[length(accepted)+1L]]<-x
    audit<-rbind(audit,data.table(File=f,Mapped=nrow(x)))
  }
  if(!length(accepted))return(list(mapping=data.table(msp_id=character(),species_name=character(),Source=character()),audit=audit))
  x<-rbindlist(accepted)
  # Conflicting non-identical names are not chosen arbitrarily.
  conflicts<-x[,.(n_names=uniqueN(species_name)),by=msp_id][n_names>1,msp_id]
  x<-x[!msp_id %in% conflicts]
  setorder(x,Source,msp_id)
  x<-unique(x,by='msp_id')
  list(mapping=x,audit=audit)
}
tax <- get_taxonomy()
taxonomy <- tax$mapping
fwrite(tax$audit,file.path(OUT,'taxonomy_source_audit.csv'))
fwrite(taxonomy,file.path(OUT,'validated_msp_taxonomy.csv'))

priority_file <- file.path(ROOT,'results','stage4','heatmap_priority_msp.csv')
feature_selection <- 'Mapped taxa chosen descriptively by detection prevalence; not a Stage 4 ranking'
priority <- character()
if(file.exists(priority_file)) {
  pri <- fread(priority_file)
  stopifnot('msp_id'%in%names(pri))
  priority <- unique(as.character(pri$msp_id))
  if(!all(priority%in%colnames(A)))stop('Unknown MSP IDs in priority file')
  feature_selection <- 'Verified saved Stage 4 priority, then mapped prevalent taxa'
}
feature_prevalence <- colMeans(A>0)
ranked <- names(sort(feature_prevalence,decreasing=TRUE))
named_ids <- taxonomy$msp_id
features <- head(unique(c(intersect(priority,named_ids),intersect(ranked,named_ids))),35L)
feature_labels <- setNames(taxonomy$species_name[match(features,taxonomy$msp_id)],features)
if(length(features)) {
  # Duplicate validated species names are made unambiguous without printing MSP IDs.
  feature_labels <- setNames(make.unique(unname(feature_labels),sep=' (variant '),features)
}
if(!length(features))message('NO VERIFIED TAXONOMY: species plots will be skipped rather than label MSP identifiers as species.')

# Create a deterministic study/group-aware representative sample set.
# Sampling does not depend on the microbial profiles.
# Allocate each group quota by largest remainder across source projects.
allocate <- function(tab,target) {
  x<-as.integer(tab); names(x)<-names(tab)
  q<-target*x/sum(x); a<-pmin(floor(q),x)
  while(sum(a)<target) {avail<-which(a<x);j<-avail[order(-(q[avail]-a[avail]),names(a)[avail])][1];a[j]<-a[j]+1}
  a
}
chosen <- character()
for(gr in c('healthy','CRC')) {
  temp<-core[group==gr]; sizes<-table(temp$study_accession)
  n<-min(150L,nrow(temp)); take<-allocate(sizes,n)
  for(st in names(take)) if(take[[st]]>0) chosen<-c(chosen,sample(temp[study_accession==st,sample],take[[st]]))
}
chosen_md <- core[match(chosen,sample)]
chosen_md <- chosen_md[order(group,study_accession,sample)]
core_ord<-core[order(group,study_accession,sample)]

to_heat <- function(samples,ids=features) {
  X <- t(relative[match(samples,rownames(relative)),ids,drop=FALSE])
  positives <- relative[,ids,drop=FALSE]
  positive_min <- min(positives[positives>0]); pseudo<-positive_min/2
  Y<-log10(X+pseudo)
  rownames(Y)<-unname(feature_labels[ids]);colnames(Y)<-samples
  Y
}
# Match scales for both heatmaps using the COMPLETE source-project-eligible population.
HEAT_BREAKS <- NULL
if(length(features)) {
  HEAT_BREAKS <- as.numeric(quantile(to_heat(core_ord$sample),c(0,.50,.98),na.rm=TRUE))
  if(length(unique(HEAT_BREAKS))<3L) { z <- range(to_heat(core_ord$sample)); HEAT_BREAKS <- c(z[1],mean(z),z[2]) }
}

# FIGURES 01-05: sample structure, alpha diversity, beta diversity.
with_fig('Figure_01','Discovery data and clinical-group composition',function(){
 x<-md0[,.(N=.N),by=.(study_accession,group)]
 p<-ggplot(x,aes(x=reorder(study_accession,N,FUN=sum),y=N,fill=group))+geom_col(width=.76)+
 geom_text(aes(label=ifelse(N>=8,N,'')),position=position_stack(vjust=.5),colour='black',fontface='bold',size=3.8)+
 coord_flip()+scale_fill_manual(values=PALETTE,drop=FALSE)+labs(title='Clinical composition across source projects',x='Source-project accession',y='Samples',fill='Clinical group')+base_theme()
 save_plot('Figure_01',p,12,8)
})
with_fig('Figure_02','Species-level richness by clinical group',function(){
 x<-copy(core);x[,richness:=rowSums(A>0)]
 p<-ggplot(x,aes(x=group,y=richness,fill=group))+geom_violin(width=.80,alpha=.30,colour=NA)+
 geom_boxplot(width=.20,outlier.alpha=.15,alpha=.88)+scale_fill_manual(values=PALETTE)+
 labs(title='Detected MSP richness across CRC and healthy samples',subtitle='Descriptive sample distribution; cohort-aware inference was evaluated separately',x=NULL,y='Detected MSP features')+base_theme()+theme(legend.position='none')
 save_plot('Figure_02',p,10,7)
})
with_fig('Figure_03','Shannon diversity by clinical group',function(){
 x<-copy(core);x[,Shannon:=vegan::diversity(A,index='shannon')]
 p<-ggplot(x,aes(x=group,y=Shannon,fill=group))+geom_violin(width=.8,alpha=.28,colour=NA)+
 geom_boxplot(width=.2,outlier.alpha=.15,alpha=.88)+scale_fill_manual(values=PALETTE)+
 labs(title='Shannon diversity across CRC and healthy samples',subtitle='Descriptive distribution; not evidence of a significant group effect',x=NULL,y='Shannon diversity')+base_theme()+theme(legend.position='none')
 save_plot('Figure_03',p,10,7)
})
with_fig('Figure_04','Bray-Curtis community ordination',function(){
 D<-vegan::vegdist(A,'bray'); ord<-stats::cmdscale(D,k=2,eig=TRUE,add=TRUE)
 x<-copy(core);x[,PCoA1:=ord$points[,1]];x[,PCoA2:=ord$points[,2]]
 p<-ggplot(x,aes(PCoA1,PCoA2,colour=group))+geom_point(size=1.3,alpha=.35)+
 scale_colour_manual(values=PALETTE)+labs(title='Bray-Curtis ordination of discovery communities',subtitle='Substantial overlap is expected; group R-squared is small',x='PCoA axis 1',y='PCoA axis 2',colour='Clinical group')+base_theme()
 save_plot('Figure_04',p,11,8)
})
with_fig('Figure_05','Source-project beta-diversity effect sizes',function(){
 f<-file.path(ROOT,'results','stage9','stage9_4_corrected_study_results.csv');if(!file.exists(f))stop('Stage 9.4 corrected result table not found')
 x<-fread(f);stopifnot(all(c('Study','R2','PERMANOVA_FDR','Dispersion_FDR')%in%names(x)))
 x[,Dispersion:=ifelse(Dispersion_FDR<.05,'Dispersion difference detected','No dispersion difference detected')]
 p<-ggplot(x,aes(x=reorder(Study,R2),y=100*R2,fill=Dispersion))+geom_col(width=.72)+coord_flip()+
 scale_fill_manual(values=c('Dispersion difference detected'='#BF6170','No dispersion difference detected'='#3D92A1'))+
 labs(title='CRC-associated Bray-Curtis variation by source project',subtitle='Source-project PERMANOVA R-squared; dispersion tests shown separately',x='Source-project accession',y='PERMANOVA R-squared (%)',fill='Dispersion test')+base_theme()
 save_plot('Figure_05',p,12,8)
})

# FIGURES 06-07: requested two versions of identical taxa-vs-sample heatmap.
with_fig('Figure_06','Prioritised microbial MSPs in 300 representative samples',function(){
 if(!length(features))stop('Validated MSP-to-species taxonomy mapping missing; provide results/stage4/msp_taxonomy.csv')
 heat_save('Figure_06',to_heat(chosen_md$sample),chosen_md,
 'Gut microbial profiles | 300 study-aware representative samples',w=18,h=10,breaks=HEAT_BREAKS)
})
with_fig('Figure_07','Prioritised microbial MSPs in all eligible discovery samples',function(){
 if(!length(features))stop('Validated MSP-to-species taxonomy mapping missing; provide results/stage4/msp_taxonomy.csv')
 heat_save('Figure_07',to_heat(core_ord$sample),core_ord,
 'Gut microbial profiles | Complete two-group discovery population',w=30,h=11,breaks=HEAT_BREAKS)
})

# FIGURES 08-10: descriptive taxa patterns, with no invented effect estimates.
with_fig('Figure_08','MSP detection frequencies by clinical group',function(){
 if(!length(features))stop('Validated MSP-to-species taxonomy mapping missing; provide results/stage4/msp_taxonomy.csv')
 ids<-head(features,20);x<-rbindlist(lapply(c('healthy','CRC'),function(g){
  ix<-core$group==g;data.table(msp_id=ids,group=g,prevalence=colMeans(A[ix,ids,drop=FALSE]>0))}))
 x[,label:=feature_labels[msp_id]]
 p<-ggplot(x,aes(x=prevalence,y=reorder(label,prevalence,FUN=mean),colour=group))+
 geom_point(size=3,position=position_dodge(width=.55))+scale_colour_manual(values=PALETTE)+
 scale_x_continuous(labels=scales::percent,limits=c(0,1))+
 labs(title='Detection prevalence of selected microbial features',subtitle=paste('Features:',feature_selection),x='Fraction of samples with nonzero signal',y=NULL,colour='Clinical group')+base_theme(12)
 save_plot('Figure_08',p,14,10)
})
with_fig('Figure_09','Within-project microbial detection frequency matrix',function(){
 if(!length(features))stop('Validated MSP-to-species taxonomy mapping missing; provide results/stage4/msp_taxonomy.csv')
 ids<-head(features,22)
 x<-rbindlist(lapply(both,function(s){
  j<-which(core$study_accession==s)
  rbindlist(lapply(c('healthy','CRC'),function(g){k<-j[core$group[j]==g];data.table(Study=s,Group=g,MSP=ids,Detection=colMeans(A[k,ids,drop=FALSE]>0))}))}))
 x[,Feature:=feature_labels[MSP]]
 p<-ggplot(x,aes(x=Study,y=Feature,fill=Detection))+geom_tile(colour='white',linewidth=.12)+
 facet_wrap(~Group,ncol=1)+scale_fill_gradient(low='#F2F6F9',high='#1B5675',limits=c(0,1),labels=percent)+
 labs(title='Microbial detection profiles by source project',subtitle='Detection frequency, not differential-abundance effect estimates',x='Source-project accession',y=NULL,fill='Detected')+base_theme(10)+
 theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
 save_plot('Figure_09',p,16,11)
})
with_fig('Figure_10','Microbial prevalence contrasts across source projects',function(){
 if(!length(features))stop('Validated MSP-to-species taxonomy mapping missing; provide results/stage4/msp_taxonomy.csv')
 ids<-head(features,20)
 x<-rbindlist(lapply(both,function(s){
  j<-which(core$study_accession==s); a<-colMeans(A[j[core$group[j]=='CRC'],ids,drop=FALSE]>0)
  b<-colMeans(A[j[core$group[j]=='healthy'],ids,drop=FALSE]>0)
  data.table(Study=s,MSP=ids,Detection_difference=100*(a-b))}))
 x[,Feature:=feature_labels[MSP]]
 lim<-max(abs(x$Detection_difference),na.rm=TRUE);lim<-max(lim,1)
 p<-ggplot(x,aes(x=Study,y=Feature,fill=Detection_difference))+geom_tile(colour='white',linewidth=.12)+
 scale_fill_gradient2(low='#176F89',mid='#F6F8F9',high='#AE3F59',midpoint=0,limits=c(-lim,lim))+
 labs(title='Source-project differences in microbial detection',subtitle='CRC minus healthy prevalence; descriptive estimates only',x='Source-project accession',y=NULL,fill='Difference (pp)')+base_theme(10)+
 theme(axis.text.x=element_text(angle=45,hjust=1))
 save_plot('Figure_10',p,16,9)
})

# FIGURES 11-12: clinically grouped relative-signal views and study sampling.
with_fig('Figure_11','Selected microbial prevalence differences',function(){
 if(!length(features))stop('Validated MSP-to-species taxonomy mapping missing; provide results/stage4/msp_taxonomy.csv')
 ids<-head(features,25); x<-data.table(MSP=ids,
  Difference=100*(colMeans(A[core$group=='CRC',ids,drop=FALSE]>0)-colMeans(A[core$group=='healthy',ids,drop=FALSE]>0)))
 x[,Feature:=feature_labels[MSP]]
 p<-ggplot(x,aes(x=reorder(Feature,Difference),y=Difference,fill=Difference>0))+
 geom_col(width=.75)+coord_flip()+scale_fill_manual(values=c('TRUE'=POS,'FALSE'=NEG),guide='none')+
 geom_hline(yintercept=0,colour='#4F5965')+
 labs(title='Descriptive CRC versus healthy detection differences',subtitle='Not a substitute for covariate-adjusted Stage 4 effect estimates',x=NULL,y='Detection prevalence difference (percentage points)')+base_theme(11)
 save_plot('Figure_11',p,13,10)
})
with_fig('Figure_12','Microbial abundance distribution among detected samples',function(){
 if(!length(features))stop('Validated MSP-to-species taxonomy mapping missing; provide results/stage4/msp_taxonomy.csv')
 ids<-head(features,8); x<-rbindlist(lapply(ids,function(id){
  y<-A[,id];data.table(MSP=id,Feature=feature_labels[id],Group=core$group,Abundance=y)}))
 x<-x[Abundance>0];x[,Log10_abundance:=log10(Abundance)]
 p<-ggplot(x,aes(x=Group,y=Log10_abundance,fill=Group))+
 geom_boxplot(outlier.alpha=.12,width=.65)+facet_wrap(~Feature,scales='free_y',ncol=4)+
 scale_fill_manual(values=PALETTE)+labs(title='Relative abundance among detected samples',subtitle='Zero-signal samples excluded; prevalence shown in other figures',x=NULL,y='Log10 relative abundance')+base_theme(11)+theme(legend.position='none')
 save_plot('Figure_12',p,15,9)
})

# Optional functional source. File is module x sample, like the MSP matrix.
ffile <- file.path(ROOT,'data','functional_signal_2340_CRC_cohort_20240329.tab')
functional <- NULL
if (file.exists(ffile)) {
 tryCatch({
  fd <- fread(ffile,showProgress=FALSE)
  sample_ids_all <- as.character(md0$sample)
  identifier <- names(fd)[1]
  # Some exports have modules as rows and samples as columns; others are transposed.
  # Identify orientation from exact sample IDs, not a guessed first-column name.
  if(all(sample_ids_all %in% names(fd))) {
    module_ids <- as.character(fd[[identifier]])
    stopifnot(!anyDuplicated(module_ids),!anyNA(module_ids))
    block <- as.matrix(fd[, ..sample_ids_all])
    storage.mode(block) <- 'double'
    F <- t(block)
    rownames(F) <- sample_ids_all
    colnames(F) <- module_ids
    functional_orientation <- 'module rows, sample columns'
  } else if(all(sample_ids_all %in% as.character(fd[[identifier]]))) {
    ord <- match(sample_ids_all,as.character(fd[[identifier]]))
    module_names <- setdiff(names(fd),identifier)
    F <- as.matrix(fd[ord,..module_names])
    storage.mode(F) <- 'double'
    rownames(F) <- sample_ids_all
    colnames(F) <- module_names
    functional_orientation <- 'sample rows, module columns'
  } else {
    stop('Functional table matches neither orientation; inspect first column and header')
  }
  stopifnot(all(is.finite(F)),all(F>=0),!anyDuplicated(colnames(F)))
  rs <- rowSums(F)
  rs[rs == 0] <- 1
  functional <- sweep(F,1,rs,'/')
 },error=function(e)message('Functional source could not be interpreted: ',conditionMessage(e)))
}
# Functional module codes are identifiers, not taxon names; display readable
# module descriptions when an authoritative mapping is available.
module_annotation <- file.path(ROOT,'results','stage5','module_annotations.csv')
module_labels <- NULL
if(file.exists(module_annotation)) {
 ma<-fread(module_annotation)
 if(all(c('module_id','module_name')%in%names(ma))) {
  ma<-unique(ma[!is.na(module_id)&!is.na(module_name)],by='module_id')
  module_labels<-setNames(as.character(ma$module_name),ma$module_id)
 }
}
# Descriptions below are from previously reported Project 4 Stage 5 outputs.
# They are NOT inferred solely from the module identifier.
verified_module_descriptions <- c(
  MGB017='Indole-3-propionic acid synthesis',
  MF0121='Propionate production via acrylate pathway',
  MF0009='Tryptophan degradation',
  MGB005='Tryptophan biosynthesis',
  MF0116='Butyrate transferase pathway',
  M00632='Galactose degradation'
)
if(is.null(module_labels)) module_labels <- verified_module_descriptions else {
  module_labels <- c(module_labels,verified_module_descriptions[!names(verified_module_descriptions) %in% names(module_labels)])
}
# Prefer the validated Stage 5 biological annotations for every available module.
stage5_annotations <- file.path(ROOT,'results','stage5','CRC_functional_modules_all.csv')
if(file.exists(stage5_annotations)) {
  sa <- fread(stage5_annotations)
  if(all(c('Module_ID','Biological_implication') %in% names(sa))) {
    sa <- unique(sa[!is.na(Module_ID) & !is.na(Biological_implication) & nzchar(Biological_implication)],by='Module_ID')
    readable <- setNames(as.character(sa$Biological_implication),as.character(sa$Module_ID))
    readable <- gsub('=>','to',readable,fixed=TRUE)
    readable <- gsub('apartate','aspartate',readable,fixed=TRUE)
    readable <- gsub('IPA synthesis','Indole-3-propionic acid synthesis',readable,fixed=TRUE)
    module_labels <- c(readable,module_labels[!names(module_labels) %in% names(readable)])
  }
}
module_priority <- names(verified_module_descriptions)
module_display <- function(ids) {
 if(is.null(module_labels))return(ids)
 z<-module_labels[ids];z[is.na(z)|!nzchar(z)]<-ids[is.na(z)|!nzchar(z)]
 z<-gsub('GDP,GTP','GDP and GTP',z,fixed=TRUE)
 z<-gsub('ADP,ATP','ADP and ATP',z,fixed=TRUE)
 z<-gsub('Citrate cycle (TCA cycle, Krebs cycle)','Citrate cycle (TCA/Krebs cycle)',z,fixed=TRUE)
 z<-gsub('Pentose phosphate pathway, non-oxidative phase, fructose 6P to ribose 5P',
   'Pentose phosphate pathway (non-oxidative)',z,fixed=TRUE)
 z<-gsub('Pentose phosphate pathway, oxidative phase, glucose 6P to ribulose 5P',
   'Pentose phosphate pathway (oxidative)',z,fixed=TRUE)
 z[grepl('^Guanine ribonucleotide biosynthesis',z)]<-'Guanine ribonucleotide biosynthesis'
 z[grepl('^propionate production',z,ignore.case=TRUE)]<-'Propionate production via acrylate pathway'
 z[grepl('^Indole-3-propionic acid synthesis',z)]<-'Indole-3-propionic acid synthesis'
 z<-sub(',.*$','',z)
 make.unique(unname(z))
}
with_fig('Figure_13','Predicted microbial function heatmap in representative samples',function(){
 if(is.null(functional))stop('Valid functional module x sample table unavailable')
 fidx<-match(core$sample,rownames(functional));prev<-colMeans(functional[fidx,,drop=FALSE]>0)
 ids<-head(unique(c(intersect(module_priority,colnames(functional)),names(sort(prev,decreasing=TRUE)))),30)
 xx<-t(functional[match(chosen_md$sample,rownames(functional)),ids,drop=FALSE])
 positive<-xx[xx>0];if(!length(positive))stop('No positive functional values')
 xx<-log10(xx+min(positive)/2)
 colnames(xx)<-chosen_md$sample
 rownames(xx)<-module_display(ids)
 heat_save('Figure_13',xx,chosen_md,'Predicted functional modules | Representative samples',w=20,h=11,legend='Log10 functional signal',annotation_legend=FALSE,legend_side='left')
})
with_fig('Figure_14','Functional potential module group contrasts',function(){
 if(is.null(functional))stop('Valid functional data unavailable')
 ix<-match(core$sample,rownames(functional));B<-functional[ix,,drop=FALSE]
 ids<-head(unique(c(intersect(module_priority,colnames(B)),names(sort(colMeans(B>0),decreasing=TRUE)))),25)
 d<-100*(colMeans(B[core$group=='CRC',ids,drop=FALSE]>0)-colMeans(B[core$group=='healthy',ids,drop=FALSE]>0))
 x<-data.table(Module=module_display(ids),Difference=d)
 p<-ggplot(x,aes(x=reorder(Module,Difference),y=Difference,fill=Difference>0))+
 geom_col()+coord_flip()+scale_fill_manual(values=c('TRUE'=POS,'FALSE'=NEG),guide='none')+
 labs(title='Functional-module detection prevalence differences',subtitle='Descriptive; pathway enrichment and inference are reported separately',x='Functional module identifier',y='CRC versus healthy detection difference (percentage points)')+base_theme(12)
 save_plot('Figure_14',p,12,9)
})

# FIGURES 15-16 use ALL THREE clinical groups; label as cross-sectional.
with_fig('Figure_15','Selected MSP detection frequency across clinical groups',function(){
 if(!length(features))stop('Validated MSP-to-species taxonomy mapping missing; provide results/stage4/msp_taxonomy.csv')
 ids<-head(features,20)
 x<-rbindlist(lapply(c('healthy','adenoma','CRC'),function(g){
  j<-which(md0$group==g);data.table(Group=g,MSP=ids,Detection=colMeans(relative[j,ids,drop=FALSE]>0))}))
 x[,Feature:=feature_labels[MSP]]
 p<-ggplot(x,aes(x=Group,y=Detection,group=Feature,colour=Feature))+
 geom_line(linewidth=.55,alpha=.8)+geom_point(size=1.8)+facet_wrap(~Feature,ncol=4)+
 scale_y_continuous(labels=percent,limits=c(0,1))+
 labs(title='Microbial detection across clinical disease groups',subtitle='Cross-sectional contrasts; lines do not represent longitudinal progression',x=NULL,y='Detection prevalence')+
 base_theme(10)+theme(legend.position='none',axis.text.x=element_text(angle=35,hjust=1))
 save_plot('Figure_15',p,17,13)
})
with_fig('Figure_16','Healthy adenoma CRC microbial prevalence matrix',function(){
 if(!length(features))stop('Validated MSP-to-species taxonomy mapping missing; provide results/stage4/msp_taxonomy.csv')
 ids<-head(features,30)
 x<-rbindlist(lapply(c('healthy','adenoma','CRC'),function(g){
  j<-which(md0$group==g);data.table(Group=g,MSP=ids,Detection=colMeans(relative[j,ids,drop=FALSE]>0))}))
 x[,Feature:=feature_labels[MSP]];x[,Group:=factor(Group,levels=c('healthy','adenoma','CRC'))]
 p<-ggplot(x,aes(x=Group,y=Feature,fill=Detection))+geom_tile(colour='white',linewidth=.7)+
 scale_fill_gradient(low='#F3F7FA',high='#1E5C76',labels=percent,limits=c(0,1))+
 labs(title='Microbial detection matrix across healthy, adenoma and CRC',subtitle='Cross-sectional detection fractions; not evidence of longitudinal progression',x='Clinical group',y=NULL,fill='Detected')+base_theme(12)
 save_plot('Figure_16',p,11,12)
})

# FIGURES 17-18 use source-project beta diversity and external validation snapshots.
with_fig('Figure_17','Discovery and external community association effect sizes',function(){
 discovery_file<-file.path(ROOT,'results','stage9','stage9_4_corrected_beta_diversity.rds')
 external_file<-file.path(ROOT,'results','stage9','beta_diversity_results.rds')
 if(!file.exists(discovery_file) || !file.exists(external_file))stop('Saved Stage 9 PERMANOVA objects unavailable')
 discovery_results<-readRDS(discovery_file);external_results<-readRDS(external_file)
 group_r2<-function(x){
  if(!'Group' %in% rownames(x) || !'R2' %in% colnames(x))stop('Group R2 absent from saved PERMANOVA result')
  as.numeric(x['Group','R2'])
 }
 x<-data.table(
  Dataset=c('Discovery (study adjusted)','Discovery (+ age, sex, BMI)',
    'External (cohort adjusted)','External (+ age, sex, BMI)'),
  R2=c(group_r2(discovery_results$primary),group_r2(discovery_results$adjusted),
    group_r2(external_results$primary),group_r2(external_results$adjusted)),
  Sample_size=c(nrow(discovery_results$metadata),nrow(discovery_results$adjusted_metadata),
    nrow(external_results$sample_metadata),nrow(external_results$adjusted_sample_metadata)),
  Stage=c('Discovery','Discovery','External validation','External validation'))
 stopifnot(all(is.finite(x$R2)),all(x$R2>0),identical(x$Sample_size,c(1868L,1637L,769L,728L)))
 p<-ggplot(x,aes(x=reorder(Dataset,R2),y=100*R2,fill=Stage))+geom_col(width=.68)+coord_flip()+
 scale_fill_manual(values=c('Discovery'='#3E8394','External validation'='#AA4D65'))+
 geom_text(aes(label=sprintf('%.3f%%',R2*100)),hjust=-.08,size=4,fontface='bold',colour='black')+
 expand_limits(y=max(x$R2*100)*1.25)+
 labs(title='Small community-composition associations in both datasets',subtitle='R-squared from completed PERMANOVA analyses; descriptive comparison, not a pooled estimate',x=NULL,y='CRC-associated variation (%)',fill='Dataset')+base_theme(12)
 save_plot('Figure_17',p,12,7)
})
with_fig('Figure_18','Discovery cohort effect and dispersion significance',function(){
 f<-file.path(ROOT,'results','stage9','stage9_4_corrected_study_results.csv');if(!file.exists(f))stop('Stage 9.4 table unavailable')
 x<-fread(f);x[,PERMANOVA:=ifelse(PERMANOVA_FDR<.05,'FDR < 0.05','Not significant')]
 x[,Dispersion:=ifelse(Dispersion_FDR<.05,'FDR < 0.05','Not significant')]
 z<-melt(x,id.vars='Study',measure.vars=c('PERMANOVA','Dispersion'),variable.name='Test',value.name='Significance')
 p<-ggplot(z,aes(x=Test,y=reorder(Study,Study),fill=Significance))+geom_tile(colour='white',linewidth=1)+
 scale_fill_manual(values=c('FDR < 0.05'='#AA3A57','Not significant'='#B5CAD1'))+
 labs(title='Project-level PERMANOVA and dispersion evidence',subtitle='Interpret community-location differences cautiously when dispersion also differs',x=NULL,y='Source-project accession',fill='Adjusted evidence')+base_theme(12)
 save_plot('Figure_18',p,10,9)
})

# FIGURES 19-20: independently tested species and modules, using the actually reported
# adjusted meta-analytic results from completed stages. These are explicit published
# project outputs, not recalculated from an incompatible abundance scale.
with_fig('Figure_19','Independent taxonomic validation evidence',function(){
 x<-data.table(Species=c('Parvimonas micra','Peptostreptococcus stomatis','Gemella morbillorum','Streptococcus salivarius'),
 OR=c(19.60,12.58,15.56,.581),Lower=c(5.45,2.11,2.47,.418),Upper=c(70.48,75.11,98.12,.806),
 FDR=c(.0306,.0306,.0306,.0306))
 p<-ggplot(x,aes(y=reorder(Species,OR),x=OR))+geom_vline(xintercept=1,linetype='dashed',colour='grey45')+
 geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.15,colour='#526977',linewidth=.8)+
 geom_point(size=3.5,colour=POS)+scale_x_log10()+
 labs(title='External validation of selected microbial detection associations',
 subtitle='Four-cohort adjusted random-effects analysis; FDR < 0.05 (reported Stage 7.9)',
 x='CRC detection odds ratio (log scale)',y=NULL)+base_theme(13)
 save_plot('Figure_19',p,12,7)
})
with_fig('Figure_20','Independent functional pathway validation evidence',function(){
 x<-data.table(Pathway=c('Propanoate fermentation via acrylate (PWY-5494)','D-galactose degradation I (PWY-6317)'),
 Estimate=c(.3942,-.1832),Lower=c(-.6349,-.7098),Upper=c(1.4232,.3435),FDR=c(.3492,.3492))
 p<-ggplot(x,aes(y=reorder(Pathway,Estimate),x=Estimate))+geom_vline(xintercept=0,linetype='dashed',colour='grey45')+
 geom_segment(aes(x=Lower,xend=Upper,yend=Pathway),linewidth=1,colour='#637784')+
 geom_point(size=4,colour='#226D89')+
 labs(title='Independent pathway testing did not establish replication',
 subtitle='Two biologically corresponding pathways; adjusted meta-analytic estimates and 95% CI',
 x='Adjusted association estimate',y=NULL)+base_theme(12)
 save_plot('Figure_20',p,12,6)
})

# FURTHER STANDALONE VISUALS: clinical provenance, geography, taxon patterns.
with_fig('Figure_21','Discovery sample counts by country',function(){
 x<-md0[,.(N=.N),by=.(country,group)]
 x[is.na(country)|!nzchar(country),country:='Not recorded']
 p<-ggplot(x,aes(x=reorder(country,N,FUN=sum),y=N,fill=group))+geom_col(width=.75)+
  geom_text(aes(label=ifelse(N>=8,N,'')),position=position_stack(vjust=.5),
    colour='black',fontface='bold',size=4)+coord_flip()+
  scale_fill_manual(values=PALETTE)+
  labs(title='Geographical diversity of discovery samples',
    subtitle='Country values are taken directly from original metadata',
    x='Country (metadata code)',y='Samples',fill='Clinical group')+base_theme(13)
 save_plot('Figure_21',p,12,8)
})
with_fig('Figure_22','Source projects and countries',function(){
 x<-md0[,.(N=.N),by=.(study_accession,country)]
 x[is.na(country)|!nzchar(country),country:='Not recorded']
 p<-ggplot(x,aes(x=country,y=study_accession,size=N,colour=N))+
  geom_point(alpha=.8)+scale_size_area(max_size=13)+
  scale_colour_gradient(low='#75B4BB',high='#99475A')+
  geom_text(aes(label=N),colour='black',fontface='bold',size=3.2,show.legend=FALSE)+
  labs(title='Source-project provenance by country',
   subtitle='One bubble per project and country combination; labels are sample counts',
   x='Country (metadata code)',y='Source-project accession',size='Samples',colour='Samples')+
  base_theme(11)+theme(axis.text.x=element_text(angle=30,hjust=1))
 save_plot('Figure_22',p,14,10)
})
with_fig('Figure_23','Microbial occurrence and abundance landscape',function(){
 if(!length(features))stop('Validated MSP-to-species taxonomy unavailable')
 ids<-head(features,35)
 x<-data.table(MSP=ids,Species=unname(feature_labels[ids]),
  Detection=100*colMeans(A[,ids,drop=FALSE]>0),
  Mean_positive=vapply(ids,function(id){z<-A[,id];if(!any(z>0))return(NA_real_);mean(z[z>0])},numeric(1)))
 x<-x[is.finite(Mean_positive)&Mean_positive>0]
 p<-ggplot(x,aes(x=Detection,y=Mean_positive,label=Species))+
  geom_point(aes(size=Detection),colour='#B34C66',alpha=.75)+
  scale_y_log10(labels=scales::label_scientific())+
  labs(title='Microbial detection and abundance are distinct properties',
    subtitle='Descriptive profiles of taxonomically verified MSPs',
    x='Samples with detected signal (%)',y='Mean relative abundance among detected samples',size='Detection (%)')+
  base_theme(12)
 save_plot('Figure_23',p,13,8)
})
with_fig('Figure_24','Country-specific microbial detection landscape',function(){
 if(!length(features))stop('Validated MSP-to-species taxonomy unavailable')
 ids<-head(features,20)
 x<-rbindlist(lapply(sort(unique(core$country)),function(ctr){
  j<-which(core$country==ctr)
  data.table(Country=ctr,MSP=ids,Detected=100*colMeans(A[j,ids,drop=FALSE]>0),N=length(j))
 }))
 x[,Species:=feature_labels[MSP]]
 p<-ggplot(x,aes(x=Country,y=Species,fill=Detected))+
  geom_tile(colour='white',linewidth=.5)+
  scale_fill_gradient(low='#EDF4F7',high='#215873')+
  labs(title='Microbial detection patterns across countries',
    subtitle='Descriptive, unadjusted detection frequency; not country-specific disease effects',
    x='Country (metadata code)',y=NULL,fill='Detected (%)')+
  base_theme(11)+theme(axis.text.x=element_text(angle=40,hjust=1))
 save_plot('Figure_24',p,15,10)
})

# Write the figure manifest after figure generation.
fwrite(manifest,file.path(OUT,'figure_manifest.csv'))
fwrite(data.table(msp_id=features,display_label=unname(feature_labels[features]),selection_method=feature_selection),
 file.path(OUT,'displayed_msp_features.csv'))
fwrite(data.table(sample=chosen_md$sample,study_accession=chosen_md$study_accession,clinical_group=as.character(chosen_md$group)),
 file.path(OUT,'representative_sample_selection.csv'))
cat('\n============================================================\n')
cat('PROJECT 4 | FIGURE ATLAS | RUN SUMMARY\n')
cat('============================================================\n')
cat('Root:',ROOT,'\nOutput:',OUT,'\n')
cat('Eligible original metadata:',nrow(md0),'\n')
cat('Source-project-adjusted CRC/healthy samples:',nrow(core),'\n')
cat('Study projects:',length(both),'\n')
cat('Representative samples:',nrow(chosen_md),'\n')
cat('Selected-feature provenance:',feature_selection,'\n')
cat('Validated species names:',length(features),'\n')
cat('Taxonomy mapping sources:',nrow(tax$audit),'\n')
cat('Functional matrix available:',!is.null(functional),'\n')
cat('Figures created:',sum(manifest$Status=='CREATED'),'of',nrow(manifest),'\n\n')
print(manifest,nrows=Inf,width=180)
cat('\nDescriptive figures do not replace the fitted statistical models.\n')
cat('Species labels come from saved taxonomy mappings; unmapped MSP plots are skipped.\n')
cat('Figures 19 and 20 summarise reported Stage 7 and Stage 8 statistics.\n')
cat('============================================================\n')
