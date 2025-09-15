# plot
#! /public/home/shilulu/anaconda3/envs/R4.2.0/bin/Rscript
setwd("F:/Github/PHD_job/2_project_hearing loss/NC_revision/badgers")
library(dplyr)
library(data.table)

my_trait <- fread("badgers_demeta_out.csv")
ukb_ano <- fread("UKB_1738trait.csv")

plt_dt <- merge(my_trait, ukb_ano, by="trait_id", all.x=TRUE)

plt_dt$shape <- ifelse(plt_dt$effect_size > 0, "positive", "negative")
my_shape <- c("positive"=24, "negative"=25)

# color for dif tissue
custom_colors <- c("Cogni/Edu" = "#FADA0C", 
                  "Diag/Symp" = "#FB9B1C", 
                  "Eye/Hearing" = "#EC2A88", 
                  "Family" = "#695092", 
                  "Lifestyle" = "#479ACD", 
                  "Medication" = "#01AEAD", 
                  "Occupation" = "#5D6967",
                  "Other" = "#c6b598")   
# set Name as factor
plt_dt$trait_id <- factor(plt_dt$trait_id, levels=unique(plt_dt$trait_id[order(plt_dt$Category)]))

# restrict yaxis at 115
plt_dt$logp = ifelse(plt_dt$pvalue == 0, 115, -log10(plt_dt$pvalue))
plt_dt$logp = ifelse(plt_dt$logp > 100, 96, plt_dt$logp)

threshold <- 0.05 / 1738
hlt_dt = plt_dt %>% 
  filter(pvalue < threshold) %>% 
  arrange(pvalue) %>% 
  head(10)
hlt_dt = hlt_dt[order(pvalue)][, .SD[1], by = trait]

hlt_list = c("Overall health rating", "Irritability", "Noisy workplace", "Miserableness", "Mood swings", "Frequency of tenseness / restlessness in last 2 weeks",
             "Frequency of depressed mood in last 2 weeks","Sleeplessness / insomnia", "Guilty feelings", "Loneliness", "Body mass index (BMI)",
             "General happiness with own health", "Health satisfaction", "Variation in diet", "Leisure/social activities: Sports club or gym")

hlt_sub = plt_dt[trait %in% hlt_list, ]
hlt_sub = hlt_sub[order(pvalue)][, .SD[1], by = trait]

hlt = rbind(hlt_dt, hlt_sub)


library(ggrepel)
library(ggplot2)
p <- ggplot(plt_dt, aes(x=trait_id, y=logp, shape=factor(shape), color=Category, fill=Category)) +
  geom_point(size=4, alpha=0.8)+
  geom_hline(yintercept=-log10(threshold), linetype="dashed", color="red", size=1)+
  scale_fill_manual(values=custom_colors) +
  scale_color_manual(values=custom_colors) +
  scale_shape_manual(values=my_shape)+
  labs(x="Traits", y=expression(-log[10](Pvalue)))+
  theme_classic()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size=14, color="black"),
        axis.title = element_text(size=14, color="black"))+
  scale_y_continuous(expand=c(0,1), limits=c(0, 100))+
  guides(fill=guide_legend(override.aes=list(size=5, shape=21)), shape="none")+
  geom_text_repel(
    data = hlt,
    aes(label = trait),
    size = 4.5,
    color="black",
    max.overlaps = Inf,
    max.iter = 5000, 
    force = 2,
    force_pull = 0.01,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50")
p

dir="F:/Github/PHD_job/2_project_hearing loss/NC_revision/badgers"
ggsave(file.path(dir,"demeta_badger_traits_add_annote.png"), p, width=15, height=6, dpi=600)
ggsave(file.path(dir,"demeta_badger_traits_add_annote.pdf"), p, width=15, height=6)
