
#----------------- plot cauchy combination p value ----//
library(ggplot2)
library(data.table)

data <- fread("gsMap_cauchy_combination_result.txt")
data <- data[order(data$P_Cauchy, decreasing=TRUE),]

# color for dif tissue
custom_colors <- c("Hair cells" = "#DD7694", 
                "Supporting cells" = "#BCD4E7", 
                "Surrounding structures" = "#056E83", 
                "Lateral wall" = "#E9D9BF", 
                "Circulating cells" = "#D4920A", 
                "Glial cells" = "#5AA4AE", 
                "Neurons" = "#65472F")   

# set Name as factor
data$Annotation <- factor(data$Annotation, levels=unique(data$Annotation))

p=ggplot(data, aes(x=Annotation, y=-log10(P_Cauchy), fill=Annotation)) +
    geom_bar(stat="identity", position="dodge", width=0.8, fill="#056E83")+
    scale_y_continuous(expand=c(0,0), limits=c(0, 9))+
    labs(x=NULL, y="-log10(P Cauchy)", title="")+
    coord_flip()+
    theme_bw()+
    theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        legend.position="none",
        legend.title=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        strip.text = element_blank(),
        strip.background=element_blank())
ggsave("plot cauchy combination p value.png", p, width=3.5, height=5, dpi=500)



#----------------- plot gsMap result for user DIY ----//
library(ggplot2)
library(RColorBrewer)
library(data.table)

data <- fread("E16.5_E1S1.MOSTA_ARHL_gsMap_plot.csv")
my_color <- rev(brewer.pal(5, "RdBu"))


"#4575b4","#74add1","#abd9e9", "#e0f3f8", "#fee090", "#fdae61", "#f46d43", "#d73027"
# plot spatial RNA alta
p = ggplot(data, aes(x = sx, y = sy)) +
    geom_point(aes(color = logp)) + 
    scale_color_gradientn(colors = c("#39489f","#39bbec","#F7F7F7", "#F4A582","#b81f25"),
                        limits = c(0, max(data$logp)),
                        breaks = seq(0, max(data$logp), by = 5)) +
    theme_minimal() +
    labs(x = "", y = "", color = expression(atop(-log[10], "(" * italic(P) * " value)"), vjust=0.2)) +
    theme(panel.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank(),
        panel.border = element_blank()) +
    theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
    theme(legend.position = c(0,0),
        legend.justification = c(0, 0),
        legend.text = element_text(color = "white", size=12),
        legend.title = element_text(color = "white", size=12),
        legend.direction = "horizontal",
        legend.key.size=unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.title.position = "top",
        legend.title.align = 0.5)

ggsave("E16.5_E1S1.MOSTA_ARHL_gsMap_plot.png", p, width=3.5, height=6, dpi=500)

# plot tissue region
my_color <-  c('Adipose tissue' = '#6b30e5', 'Adrenal gland'='#bf024fff', 'Bone'='#103a14ff', 'Brain'='#ef833aff', 
            'Cartilage'='#35586dff', 'Cartilage primordium'='#3cb44bff', 'Cavity'='#dfdce0ff', 
            'Choroid plexus'='#bd3addff', 'Connective tissue'='#0bd3b1ff', 'Dorsal root ganglion'='#b74c11ff', 
            'Epidermis'='#036df4ff', 'GI tract'='#5c5ca6ff', 'Heart'='#d3245aff', 'Inner ear'='#03fff4ff', 
            'Jaw and tooth'='#f062f9ff', 'Kidney'='#62cfe8ff', 'Liver'='#c923b1ff', 'Lung'='#7ec136ff', 
            'Meninges'='#dfca43ff', 'Mucosal epithelium'='#2f7dd1ff', 'Muscle'='#af1041ff', 'Smooth muscle'='#fc5151ff', 
            'Spinal cord'='#f9d5baff', 'Submandibular gland'='#ab32e6ff', 'Sympathetic nerve'='#cc5a0dff')

p = ggplot(data, aes(x = sx, y = sy)) +
    geom_point(aes(color = annotation)) + 
    scale_color_manual(values=my_color) +
    theme_minimal() +
    labs(title = "", x = "", y = "") +
    theme(panel.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank(),
        panel.border = element_blank()) +
    theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
    theme(legend.position = "right",
        legend.text = element_text(color = "black"),
        legend.title = element_blank()) 
ggsave("E16.5_E1S1.MOSTA_ARHL_cell_region.png", p, width=7, height=6, dpi=500)