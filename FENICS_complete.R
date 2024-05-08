library(librarian)
shelf(readxl, Hmisc, tidyverse, ontologyIndex, magrittr, ggrepel, quiet=T)

# Tree Setup ----
## OBO to Is_A and Ancs ----
raw_tree <- get_OBO("fenics_v8.obo", extract_tags="everything")
fenics_ids <- tibble(raw_id = raw_tree$id,
                     fenics_id = sapply(raw_tree$property_value,
                                        function(x){
                                          str_extract(paste(x,sep="",collapse=";"),"FENICS:\\d+")}))
for (i in c("id", "parents", "children", "ancestors")) {
  for (j in 1:length(raw_tree[[i]])) {
    for (k in 1:length(raw_tree[[i]][[j]])) {
      l <- which(fenics_ids$raw_id==raw_tree[[i]][[j]][k])
      raw_tree[[i]][[j]][k] <- fenics_ids$fenics_id[l]
    }
    raw_tree[[i]][[j]] <- paste(raw_tree[[i]][[j]], sep="",collapse=";")
  }
  raw_tree[[i]] <- as.character(raw_tree[[i]])
}
raw_tibble <- tibble(term.id=raw_tree$id, term=raw_tree$name, parents=raw_tree$parents,
                     ancestors=raw_tree$ancestors, children=raw_tree$children) %>%
  filter(term.id != "NA")
isa <- raw_tibble %>% select(term.id, term, is.a=parents) %>%
  separate_rows(is.a, sep=";")
def <- raw_tibble %>% select(term.id, term) %>% unique()
ancs <- raw_tibble %>% select(term.id, term, ancestors) %>%
  mutate(ancestors=paste(term.id, ancestors, sep=";")) %>%
  separate_rows(ancestors, sep=";") %>% unique
conv_dict <- fenics_ids %>% mutate(raw_id = gsub("http://webprotege.stanford.edu/", "", raw_id))
comps <- ancs %>%
  filter(ancestors %in% c("FENICS:0004", "FENICS:0003",
                          "FENICS:0006")) %>%
  mutate(effect = case_when(ancestors=="FENICS:0004" ~ "lof",
                            ancestors=="FENICS:0003" ~ "gof",
                            ancestors=="FENICS:0006" ~ "mix")) %>%
  select(-ancestors) %>%
  unique()
cats <- paste0("FENICS:", c("0009","0016","0021","0033","0038","0045","0052","0059","0071","0075",
                            "0078","0079","0082","0097","0102","0106","0110","0115","0119","0123",
                            "0128","0132"))
param <- ancs %>%
  filter(ancestors %in% cats) %>%
  unique()
params <- def %>%
  filter(term.id %nin% param$term.id) %>%
  mutate(ancestors="NA") %>%
  bind_rows(param) %>%
  rename(def=term) %>%
  left_join(def, by=c("ancestors"="term.id")) %>%
  mutate(parameter=case_when(is.na(term) ~ "other", T ~ gsub("Effect on ", "", term))) %>%
  select(term.id, term=def, parameter)

# Analysis ----
## Descriptive Stats ----
###PER Dictionary = paralogous mapping for SCN1A/2A/3A/8A ... 
###... obtained from Perez-Palma et al., 2020 (31871067) ...
###... expanded to 3 columns: index, gene, position 
per <- read_csv("per_dict.csv")
raw <- read_csv("Parthasarathy2024_AJHG_Table_S1.csv") %>% 
  filter(grepl("SCN", gene)) %>% 
  unite("expID", c(gene, variant, PMID), remove = FALSE)
prop <- raw %>% filter(!grepl("Normal",term), !grepl(" normal", term)) %>%
  left_join(ancs) %>% select(-term, -term.id) %>%
  rename(term.id=ancestors) %>% left_join(def) %>%
  bind_rows(raw %>% filter(grepl("Normal",term) | grepl(" normal",term))) %>%
  left_join(comps) %>% distinct()
overall_fxn <- prop %>% filter(term.id %in% c("FENICS:0137","FENICS:0145","FENICS:0141")) %>%
  mutate(overall = case_when(term.id=="FENICS:0137" ~ "GoF",
                             term.id=="FENICS:0145" ~ "Mix",
                             term.id=="FENICS:0141" ~ "LoF")) %>%
  select(expID, overall) %>% distinct()
prop %<>% left_join(overall_fxn)
raw %>% pull(expID) %>% unique %>% length
raw %>% select(gene, expID) %>% unique %>% count(gene)
raw %>% select(gene, variant) %>% unique %>% nrow
raw %>% select(gene, variant) %>% unique %>% count(gene)
raw %>% nrow
prop %>% select(expID, overall) %>% unique() %>% count(overall)
n_measurable <- raw %>% filter(expID %nin% (raw %>% filter(grepl("FENICS:0083",term.id)) %>% pull(expID))) %>% pull(expID) %>% unique %>% length
opp_features <- prop %>% left_join(comps) %>% filter(!is.na(effect), effect != "mix") %>% select(expID,effect) %>% unique %>% count(expID) %>% filter(n>1)
nrow(opp_features)
nrow(opp_features)/n_measurable
opp_features %>% left_join(overall_fxn) %>% count(overall)
overall_fxn %>% count(overall)
prop %>% filter(overall=="GoF") %>% unique() %>%
  count(term) %>% arrange(-n)
prop %>% filter(overall=="Mix") %>% unique() %>%
  count(term) %>% arrange(-n) %>% View
prop %>% filter(overall=="LoF") %>% unique() %>%
  count(term) %>% arrange(-n) %>% View
raw %>% count(expID) %>% pull(n) %>% range
raw %>% count(expID) %>% pull(n) %>% median
raw %>% pull(term.id) %>% unique %>% length
pos <- raw %>% filter(!grepl("Normal",term), !grepl(" normal", term))
pos %>% nrow
pos %>% pull(term.id) %>% unique %>% length
pos %>% count(expID) %>% pull(n) %>% median
pos %>% count(expID) %>% pull(n) %>% range
neg <- raw %>% filter(term.id %nin% pos$term.id)
neg %>% nrow
neg %>% pull(term.id) %>% unique %>% length
neg %>% count(expID) %>% pull(n) %>% median
neg %>% count(expID) %>% pull(n) %>% range
raw %>% count(term) %>% arrange(-n) %>% mutate(n=n/length(unique(raw$expID))) %>% View
prop %>% pull(term.id) %>% unique %>% length
prop %>% count(expID) %>% pull(n) %>% median
prop %>% count(expID) %>% pull(n) %>% range
prop %>% count(term) %>% arrange(-n) %>% mutate(n=n/length(unique(raw$expID))) %>% View
ps3_cal_terms <- c("Effect on persistent current",
                   "Effect on shift of voltage dependence of fast inactivation",
                   "Effect on shift of voltage dependence of activation",
                   "Effect on peak current")
ps3_app <- prop %>% filter(term %in% ps3_cal_terms) %>% unique %>%
  select(gene,variant) %>% distinct()
raw %>% filter(variant %nin% ps3_app$variant) %>%
  count(term) %>% arrange(-n) 
prop %>% filter(variant %nin% ps3_app$variant) %>%
  filter(!grepl("Normal",term), !grepl(" normal",term)) %>%
  count(term) %>% arrange(-n) 

###Structural features of sodium channels obtained from UniProt ... 
###... with columns domain (I-IV), segment (1-6), name (e.g., I1 or III5), 
###... start and end based on SCN1A, to be then mapped to paralogous indices
regions <- read_csv("fenics_nav_topo.csv")
lollin <- prop %>% select(expID, gene, variant) %>%
  mutate(position = as.numeric(str_extract(variant, "\\d+"))) %>%
  left_join(per) %>%
  unique() %>%
  select(-expID) %>%
  group_by(gene) %>%
  count(index)
plot <- ggplot(lollin) +
  geom_segment(alpha=1, size=1, aes(x=index,xend=index,y=0,yend=n+2)) +
  geom_point(size=3,  alpha=1, color="steelblue", aes(x=index,y=n+2)) +
  geom_rect(aes(xmin=1,xmax=2188,ymin=-2,ymax=0)) +
  geom_rect(data=regions, aes(xmin=start,xmax=end,ymin=-2.2,ymax=.2,fill=domains)) +
  theme_classic() +
  scale_y_continuous(limits=c(-2.2,10)) +
  labs(x="Position Index") +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.title = element_blank()) +
  facet_grid(rows="gene", labeller=label_wrap_gen(width=10))
plot

## Heatmap and radar ----
x <- prop %>%
  expand(expID,effect) %>%
  left_join(prop) %>%
  left_join(params) %>%
  filter(parameter != "other") %>%
  filter(expID %in% (prop %>% filter(!is.na(effect)) %>% pull(expID))) %>%
  filter(!is.na(effect))
y <- prop %>%
  mutate(effect = case_when(overall=="Mix" ~ "mix",
                            overall=="GoF" ~ "gof",
                            overall=="LoF" ~ "lof")) %>%
  select(expID, effect, overall) %>% unique() %>%
  mutate(parameter = "overall effect",
         effect=factor(effect, levels=c("lof","mix","gof"))) %>%
  arrange(effect)
z <- x %>%
  bind_rows(y) %>%
  filter(!is.na(overall)) %>%
  mutate(expID = factor(expID, levels=y$expID),
         parameter = factor(parameter, levels = c("overall effect", x%>%count(parameter)%>%arrange(-n)%>%pull(parameter))))
x%>%count(parameter)%>%arrange(-n)%>%pull(parameter) %>% setNames(NULL)
ggplot(z,aes(x=parameter,y=expID,fill=effect)) + geom_tile() +
  scale_fill_manual(values=c("gof"="steelblue","lof"="tomato","mix"="#A2727D"),na.value="gray90") +
  theme_classic() +
  theme(legend.position="none", axis.line=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),axis.ticks=element_blank())
radata <- prop %>%
  filter(!grepl("Normal", term), !grepl(" normal", term)) %>%
  left_join(isa) %>% left_join(def %>% rename(is.a=term.id, is.a_def=term)) %>%
  filter(grepl("Component", is.a_def)) %>%
  group_by(overall) %>% count(term) %>% arrange(desc(n)) %>%
  filter(n>4) %>%
  ungroup()
radata %<>% expand(overall,term) %>%
  left_join(radata) %>%
  replace_na(list(n=0)) %>%
  left_join(prop %>% select(expID,overall) %>% unique %>% count(overall,name="tot")) %>%
  mutate(frac=n/tot) %>%
  left_join(comps) %>%
  group_by(effect) %>%
  arrange(desc(n), .by_group=T) %>%
  ungroup() %>%
  mutate(term = factor(term, levels=unique(term)))
lapply(c("GoF","LoF","Mix"), function(i) {
  radata %>% filter(overall==i) %>%
    ggplot() + coord_polar() +
    geom_col(aes(x=term, y=frac)) +
    geom_segment(aes(x=0, xend=22.5, y=0.75, yend=0.75), color="gray", alpha=0.5) +
    geom_segment(aes(x=0, xend=22.5, y=0.5, yend=0.5), color="gray", alpha=0.5) +
    geom_segment(aes(x=0, xend=22.5, y=0.25, yend=0.25), color="gray", alpha=0.5) +
    geom_col(aes(x=term, y=frac, fill=effect),width=1) +
    scale_y_continuous(limits=c(0,.75), expand=c(0,0)) +
    scale_fill_manual(values=c("steelblue","tomato")) +
    theme_classic() +
    theme(text=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(),
          legend.position="none",
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA))
})

## ClinGen Calibration ----
###Nested functions
###Requires table with column names: gene, variant, variant_class, ...
###... parameter, mean_adj ...
###... where variant_class is 'benign' or 'pathogenic' and ...
###... mean_adj is mean deviation from WT, adjusted to be compatible ...
###... i.e., ensuring peak current is scaled so that 0 = absent and 1 = WT ...
###... and that V1/2 values are taken as absolute value of deviation. 
###Output: list where entry 1 is final table of thresholds, and ... 
###... entry 2 records the computations in detail. 
comp_lrs <- function(ctrl) {
  func_lrs_maxs <- tibble(variant_class=c("benign","pathogenic"), dir="more",threshold=-1000)
  func_lrs_result_new <- tibble()
  func_lrs_result <- ""
  func_lrs_result_pre <- tibble(x="")
  
  while(!isTRUE(all.equal(func_lrs_result_new, func_lrs_result_pre))) {
    func_lrs_result_pre <- func_lrs_result
    func_lrs_result <- func_lrs_result_new
    ctrl_adj1 <- ctrl %>% left_join(func_lrs_maxs) %>%
      mutate(norm = case_when((dir == "less" & mean_adj < threshold) ~ "abn",
                              (dir == "more" & mean_adj >= threshold) ~ "abn",
                              (dir == "less" & mean_adj >= threshold) ~ "nrm",
                              (dir == "more" & mean_adj < threshold) ~ "nrm"),
             norm = factor(norm, levels=c("nrm","abn")))
    
    func_lrs <- lapply(unique(ctrl$parameter), lr_param, ctrl, ctrl_adj1) %>% reduce(bind_rows) %>%
      mutate(strength = case_when(LR >= 350^(1/2) ~ "Strong",
                                  LR >= 350^(1/4) ~ "Moderate",
                                  LR >= 350^(1/8) ~ "Supporting",
                                  T ~ "None"),
             strength=factor(strength, levels=c("None","Supporting","Moderate","Strong")),
             threshold = case_when(dir=="less" ~ -threshold, dir=="more" ~ threshold)) %>%
      group_by(parameter,strength, dir)
    func_lrs_result_new <- func_lrs %>%
      slice_min(threshold) %>%
      ungroup() %>%
      filter(strength!="None") %>%
      mutate(threshold=abs(threshold))
    func_lrs_maxs <- func_lrs %>%
      slice_max(LR) %>%
      ungroup() %>% group_by(parameter,dir) %>%
      slice_max(strength) %>%
      ungroup() %>%
      filter(strength!="None") %>%
      mutate(threshold=abs(threshold))
  }
  return(list(func_lrs_result_new, func_lrs))
}

lr_param <- function(p, ctrl, ctrl_adj1) {
  print(p)
  ctrl_adj <- ctrl_adj1 %>% group_by(mean_adj) %>% slice_max(norm) %>%
    filter(norm=="nrm",parameter==p) %>%
    inner_join(ctrl_adj1 %>%
                 filter(norm=="abn", parameter!=p) %>% select(variant_class,gene,variant) %>% unique)
  ctrl_param <- ctrl %>% filter(parameter==p) %>% anti_join(ctrl_adj)
  vals <- lapply(unique(ctrl_param$mean_adj), lr_from_thres, dir="more", ctrl_param=ctrl_param) %>%
    reduce(bind_rows) %>%
    bind_rows(lapply(unique(ctrl_param$mean_adj), lr_from_thres, dir="less", ctrl_param=ctrl_param) %>%
                reduce(bind_rows)) %>%
    add_column(parameter=p)
  return(vals)
}

lr_from_thres <- function(thres, dir, ctrl_param) {
  tp <- ctrl_param %>% filter(variant_class=="pathogenic") %>% filter(mean_adj < thres) %>% nrow
  fn <- ctrl_param %>% filter(variant_class=="pathogenic") %>% filter(mean_adj >= thres) %>% nrow
  tn <- ctrl_param %>% filter(variant_class=="benign") %>% filter(mean_adj >= thres) %>% nrow
  fp <- ctrl_param %>% filter(variant_class=="benign") %>% filter(mean_adj < thres) %>% nrow
  if (dir=="more") {
    swap <- tp
    tp <- fn
    fn <- swap
    swap <- tn
    tn <- fp
    fp <- swap
  }
  if(fp==0) {fp <- 1}
  if(tn==0) {tn <- 1}
  lrp1 <- (tp/(tp+fn))/(1-(tn/(tn+fp)))
  #lrp1 <- (tp*tn)/(fn*fp)
  #lrp1 <- tp/(tp+fn)
  #lrp1 <- fisher.test(matrix(c(tp,fn,fp,tn),nrow=2))$p.value
  return(tibble(threshold=thres, dir=dir, TP=tp, FN=fn,
                FP=fp, TN=tn, LR=lrp1))
}

###scn_raw table can be obtained by combining ...
###... pathogenic control data from from our Table S2 ...
###... with as yet unpublished benign control data. 
scn_raw <- read_csv("scn_fxn_ctrl_all.csv")
scn_thresh <- comp_lrs(scn_raw) 

###kcn_raw table modified from published data 
###Supplemental Tables 4 and 5 from Vanoye et al., 2022 (35104249)
###We use an 'assay' column to distinguish 'het', i.e., heterozygous ...
###... from 'hom', i.e., homozygous experimental conditions
###A quick lolliplot precedes calibration as above
kcn_raw <- read_csv("kcnq2_fxn_ctrl_all.csv")
kcn_regions <- tribble(
  ~segments, ~start, ~end,
  "1", 92, 112,
  "2", 123, 143,
  "3", 167, 187,
  "4", 196, 218,
  "5", 232, 252,
  "5-6", 264, 286,
  "6", 291, 312
)
kcn_lollin <- kcn_raw %>% select(gene, variant, variant_class) %>%
  mutate(position = as.numeric(str_extract(variant, "\\d+"))) %>%
  unique() %>%
  group_by(variant_class) %>%
  count(position)
plot <- ggplot(kcn_lollin) +
  geom_segment(alpha=1, size=1, aes(x=position,xend=position,y=0,yend=n)) +
  geom_point(size=3,  alpha=1, aes(x=position,y=n, color=variant_class)) +
  geom_rect(aes(xmin=1,xmax=872,ymin=-1,ymax=0)) +
  geom_rect(data=kcn_regions, aes(xmin=start,xmax=end,ymin=-1.1,ymax=.1,fill=segments)) +
  theme_classic() +
  scale_y_continuous(limits=c(-1.1,4.5), expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=c("pathogenic"="red","benign"="steelblue")) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.title = element_blank())
plot
kcn_het_thresh <- comp_lrs(kcn_raw %>% filter(assay == "het"))
kcn_hom_thresh <- comp_lrs(kcn_raw %>% filter(assay == "hom"))

###Plots used in Figures 3 and 4 are built using these frameworks
kcn_hom_thresh %>%
  filter(param=="Current Density", dir=="less", threshold>=0) %>%
  ggplot(aes(x=100*threshold,y=LR)) +
  annotate(geom="rect",xmin=0,xmax=100,ymin=0,ymax=350^(1/8), fill="gray80",alpha=.5) +
  annotate(geom="rect",xmin=0,xmax=100,ymin=350^(1/8),ymax=350^(1/4), fill="steelblue4",alpha=.5) +
  annotate(geom="rect",xmin=0,xmax=100,ymin=350^(1/4),ymax=350^(1/2), fill="steelblue3",alpha=.5) +
  annotate(geom="rect",xmin=0,xmax=100,ymin=350^(1/2),ymax=25, fill="steelblue1",alpha=.5) +
  geom_line(color="black",size=1.5) +
  labs(x="Percent Reduction in Peak Current", y="Likelihood Ratio") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  theme_classic()

###All dots visualizations are built by iterating the contents ...
###... of the comp_lrs and lr_param function to obtain ctrl_param ...
###... which is a table of control values for a single parameter ...
###... for a single iteration of the analysis. 
ctrl_param %>% 
  ggplot(aes(x=0, y=mean_adj)) +
  geom_rect(xmin=-.75, xmax=.75, ymin=0, ymax = 27, fill="white", color="black") +
  geom_hline(yintercept=3.2, linetype="dashed", color="red") +
  geom_hline(yintercept=2.95, linetype="dashed", color="salmon") +
  geom_jitter(aes(color=type), width=.2, size=4) +
  scale_x_continuous(limits=c(-.8,.8), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,27), expand=c(0,0)) +
  scale_color_manual(values=c("pathogenic"="red","benign"="steelblue")) +
  theme_classic() +
  theme(text=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(),
        legend.position="none")

###Benign densities in Figure 4 are built with this frame
hom_begn_v12a <- kcn_raw %>% filter(assay == "hom") %>% 
  filter(parameter=="V1/2 Activation", variant_class=="benign") %>%
  arrange(-mean_adj) %>%
  mutate(pct = row_number()/24)
ggplot(hom_begn_v12a, aes(y=mean_adj)) + geom_density(fill="steelblue",alpha=0.4) +
  coord_cartesian(ylim=c(0,27), expand=c(0,0)) +
  geom_hline(yintercept=7.2, linetype="dashed", color="red") +
  geom_hline(yintercept=5.7, linetype="dashed", color="salmon") +
  theme_classic() +
  theme(axis.title=element_blank(),axis.ticks=element_blank(),
        axis.text=element_blank(), axis.line=element_blank())


## Structure-Function Associations ----
anemonemone <- function(groupcol, intable) {
  groupnames <- intable %>% pull(!!sym(groupcol)) %>% unique() %>% na.omit()
  return(lapply(groupnames, anemone, groupcol=groupcol, intable=intable) %>% reduce(bind_rows))
}
anemone <- function(groupname, groupcol, intable) {
  groupids <- intable %>% filter(!!sym(groupcol)==groupname) %>% pull(expID)
  fxs <- intable %>% filter(expID %in% groupids) %>% pull(term.id) %>% unique
  return(lapply(fxs, nemo, groupids=groupids, groupname=groupname, intable=intable) %>%
           reduce(bind_rows) %>% left_join(ancs, by=c("term.id")) %>% select(-ancestors) %>% distinct())
}
nemo <- function(fx, groupids, groupname, intable) {
  yy <- prop %>% filter(term.id==fx, expID %in% groupids, expID %in% intable$expID) %>% pull(expID) %>% unique()
  nyy <- length(yy)
  ny <- prop %>% filter(term.id==fx, expID %nin% groupids, expID %in% intable$expID) %>% pull(expID) %>% unique()
  nny <- length(ny)
  yn <- prop %>% filter(expID %in% groupids, expID %nin% yy) %>% pull(expID) %>% unique()
  nyn <- length(yn)
  nn <- prop %>% filter(expID %nin% groupids, expID %nin% ny, expID %in% intable$expID) %>% pull(expID) %>% unique()
  nnn <- length(nn)
  ft <- fisher.test(matrix(c(nyy,nny,nyn,nnn),nrow=2))
  p <- ft$p.value
  OR <- ft$estimate
  OR_lo <- ft$conf.int[1]
  OR_hi <- ft$conf.int[2]
  if(nny==0 | nyn==0) {
    OR_adj <- (nyy+.5)*(nnn+.5)/((nny+.5)*(nyn+.5))
    se <- sqrt(1/nyy+1/nny+1/nyn+1/nnn)
    OR_adj_lo <- exp(log(OR_adj) - 1.96*se)
    OR_adj_hi <- exp(log(OR_adj) + 1.96*se)
  } else {OR_adj <- OR; OR_adj_lo <- OR_lo; OR_adj_hi <- OR_hi}
  result <- tribble(~term.id, ~group, ~n_yy, ~freq_in, ~freq_out, ~pval, ~log10p, ~OR,
                    ~OR_lo, ~OR_hi, ~OR_adj, ~OR_adj_lo, ~OR_adj_hi,
                    fx, groupname, nyy, nyy/(nyy+nyn), nny/(nny+nnn), p, -log10(p), OR,
                    OR_lo, OR_hi, OR_adj, OR_adj_lo, OR_adj_hi)
  return(result)
}
trop <- prop %>%
  filter(!grepl("Mild", term), !grepl("Moderate", term), !grepl("Severe", term)) %>%
  mutate(position=as.numeric(str_extract(variant, "\\d+"))) %>%
  left_join(read_csv("./Documents/per_dict.csv"), by=c("gene","position")) %>%
  mutate(domain = case_when(index %in% 0:138 ~ "N",
                            index %in% 138:460 ~ "I",
                            index %in% 461:828 ~ "I-II",
                            index %in% 829:1063 ~ "II",
                            index %in% 1064:1364 ~ "II-III",
                            index %in% 1365:1627 ~ "III",
                            index %in% 1628:1690 ~ "III-IV",
                            index %in% 1691:1935 ~ "IV",
                            index %in% 1936:2188 ~ "C"),
         segment = case_when(index %in% c(138:160,829:847,1365:1382,1691:1708) ~ "S1",
                             index %in% c(168:188,859:878,1396:1414,1751:1768) ~ "S2",
                             index %in% c(203:220,893:914,1429:1447,1751:1768) ~ "S3",
                             index %in% c(228:244,917:934,1456:1474,1784:1800) ~ "S4",
                             index %in% c(264:283,951:969,1492:1511,1820:1837) ~ "S5",
                             index %in% c(408:432,1007:1027,1567:1588,1860:1882) ~ "S5-6",
                             index %in% c(440:460,1043:1063,1606:1627,1913:1935) ~ "S6"))
seg_assoc <- anemonemone("segment", trop) %>%
  arrange(pval) %>%
  filter(OR_adj>1) %>%
  filter((pval!=lead(pval)) | (group!=lead(group)), !grepl("Normal", term), !grepl(" normal", term)) %>% #,
  mutate(n = row_number(), t = max(n), pcut = .1*n/t,
         FDR = case_when(pval < pcut ~ "SIG", T ~ "NOT")) %>%
  mutate(FDR = case_when(row_number() < max(which(FDR=="SIG")) ~ "SIG", T ~ FDR))
dom_assoc <- anemonemone("domain", trop) %>%
  arrange(pval) %>%
  filter(OR_adj>1) %>%
  filter((pval!=lead(pval)) | (group!=lead(group)), !grepl("Normal", term), !grepl(" normal", term)) %>% #,
  mutate(n = row_number(), t = max(n), pcut = .1*n/t,
         FDR = case_when(pval < pcut ~ "SIG", T ~ "NOT")) %>%
  mutate(FDR = case_when(row_number() < max(which(FDR=="SIG")) ~ "SIG", T ~ FDR))
seg_assoc_fig <- seg_assoc %>%
  left_join(tibble(
    group = c("S1","S2","S3","S4","S5","S5-6","S6"),
    start = c(829,859,893,917,951,1007,1043),
    end = c(847,878,914,934,969,1027,1063)
  )) %>%
  mutate(nomsig = case_when(FDR=="SIG" ~ "sig_fdr", pval < 0.05 ~ "sig_nom", T ~ "sig_not"),
         term.lab = case_when(nomsig == "sig_not" ~ "", T ~ term),
         term.lab = gsub("voltage dependence of", "V1/2",
                         gsub("Hyperpolarizing", "Left",
                              gsub("Depolarizing", "Right", term.lab)))) %>%
  arrange(-pval)
dom_assoc_fig <- dom_assoc %>%
  left_join(tibble(
    group = c("N","I","I-II","II","II-III","III","III-IV","IV","C"),
    start = c(1,138,461,829,1064,1365,1628, 1691, 1936),
    end = c(137,460,828,1063,1364,1627,1690, 1935, 2188)
  )) %>%
  mutate(nomsig = case_when(FDR=="SIG" ~ "sig_fdr", pval < 0.05 ~ "sig_nom", T ~ "sig_not"),
         term.lab = case_when(nomsig == "sig_not" ~ "", T ~ term),
         term.lab = gsub("voltage dependence of", "V1/2",
                         gsub("Hyperpolarizing", "Left",
                              gsub("Depolarizing", "Right", term.lab))),
         start_tm = case_when(group %in% c("I","II","III","IV") ~ start, T ~ 0),
         end_tm = case_when(group %in% c("I","II","III","IV") ~ end, T ~ 0)) %>%
  arrange(-pval)
ggplot(data=dom_assoc_fig, aes(x=.5*(start+end), y=log10(OR_adj))) +
  geom_rect(aes(xmin=floor(min(start)/25)*25,xmax=ceiling(max(end)/25)*25,ymin=-.12,ymax=-.04), fill="gray50") +
  geom_rect(aes(xmin=start_tm,xmax=end_tm,ymin=-.16,ymax=0, fill=group)) +
  geom_vline(aes(xintercept=end), color="gray50",linetype="dashed") +
  geom_vline(aes(xintercept=start), color="gray50",linetype="dashed") +
  geom_jitter(aes(color=nomsig, size=log10p), width=40) +
  scale_color_manual(values=c("sig_fdr"="red","sig_nom"="steelblue", "sig_not"="gray50")) +
  scale_fill_manual(values=c("I"="#F8766D","II"="#7CAE00", "III"="#00BFC4", "IV"="#C77CFF")) +
  scale_size_area(max_size=12) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(limits=c(-.16,1.5), expand=c(0,0)) +
  theme_classic() +
  theme(
    text=element_blank(),
    axis.line=element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")
ggplot(data=seg_assoc_fig, aes(x=.5*(start+end), y=log10(OR_adj))) +
  geom_rect(aes(xmin=floor(min(start)/25)*25,xmax=ceiling(max(end)/25)*25,ymin=-.09,ymax=-.03), fill="#B8D56D") +
  geom_rect(aes(xmin=start,xmax=end,ymin=-.12,ymax=0), fill="#7CAE00") +
  geom_vline(aes(xintercept=end), color="gray50",linetype="dashed") +
  geom_vline(aes(xintercept=start), color="gray50",linetype="dashed") +
  geom_jitter(aes(color=nomsig, size=log10p), width=5) +
  scale_color_manual(values=c("sig_fdr"="red","sig_nom"="steelblue", "sig_not"="gray50")) +
  scale_size_area(max_size=12) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(limits=c(-.12,1.54), expand=c(0,0), breaks=c(0,log10(3),1,log10(30))) +
  theme_classic() +
  theme(
    text=element_blank(),
    axis.line=element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")

## SCN2A Projection ----
phen <- read_csv("./Documents/SCN2A_full_V16.csv")
raw %>% filter(gene=="SCN2A") %>%
  pull(variant) %>% unique() %>% length
twoa <- raw %>% filter(gene=="SCN2A") %>% pull(variant)
phen %>% filter(variant %in% twoa) %>%
  select(-HPO, -definition, -hpo_term) %>% unique() %>% nrow
phen %>% filter(variant_type_2=="missense") %>%
  select(-HPO, -definition, -hpo_term) %>% unique() %>% nrow
per2a <- per %>% filter(gene=="SCN2A")
recur <- phen %>% filter(variant_type_2=="missense") %>% select(famID, variant) %>%
  unique() %>% count(variant) %>% filter(n>2)
sum(recur$n)
phen %>% select(famID, variant) %>% unique() %>%
  filter(variant %in% recur$variant) %>%
  filter(variant %in% twoa) %>%
  nrow
pull(variant) %>% unique %>% length
phenmis %>%
  mutate(position = as.numeric(str_extract(variant, "\\d+"))) %>%
  left_join(per2a) %>%
  filter(index %in% (raw %>%
                       mutate(position=as.numeric(str_extract(variant,"\\d+"))) %>%
                       left_join(per) %>% pull(index))) %>%
  nrow
