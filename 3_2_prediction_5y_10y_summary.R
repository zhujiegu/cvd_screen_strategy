library(ggplot2)
load(file = paste0("~/epi_paper/risk_ratio.RData"))
# save(risk_ratio, file = paste0("~/epi_paper/risk_ratio.RData"))

risk_gr_colors <- c("red4", "red","darkorange", "#ffd200", "springgreen3")
names(risk_gr_colors)<- risk_ratio$risk_class %>% levels
risk_labels <- c('>10%','7.5%-10%','5%-7.5%', '2.5%-5%','<2.5%')


p_female_ratio <- ggplot(risk_ratio %>% filter(gender=='female'), aes(x = factor(lm_age), y = risk_ratio, color=risk_class, fill = risk_class)) +
  geom_bar(alpha=0.6,stat = "identity", position = position_dodge()) +
  labs(
    x = "Landmark age",
    y = "Risk Ratio",
    title = "Women",
    color = "Risk group",  # Legend title for color
    fill = "Risk group"    # Legend title for fill
  ) +
  theme_minimal() +
  scale_color_manual(values=risk_gr_colors, labels=risk_labels) +
  scale_fill_manual(values=risk_gr_colors, labels=risk_labels) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_male_ratio <-ggplot(risk_ratio %>% filter(gender=='male'), aes(x = factor(lm_age), y = risk_ratio, color=risk_class, fill = risk_class)) +
  geom_bar(alpha=0.6,stat = "identity", position = position_dodge()) +
  labs(
    x = "Landmark age",
    y = "Risk Ratio",
    title = "Men",    
    color = "Risk group",  # Legend title for color
    fill = "Risk group"    # Legend title for fill
  ) +
  theme_minimal() +
  scale_color_manual(values=risk_gr_colors, labels=risk_labels) +
  scale_fill_manual(values=risk_gr_colors, labels=risk_labels) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p_male_ratio)

jpeg(file = '~/epi_paper/outp_figs/risk_ratio.jpeg', units="in", width=8, height=5, res=600)
grid.arrange(arrangeGrob(p_female_ratio + theme(legend.position="none"),
                         p_male_ratio + theme(legend.position="none"), nrow=2),
             mylegend, ncol=2,widths=c(10, 1.5))
dev.off()
