library(tidyverse)
library(ggridges)

data = read.csv("SX_ridgeline_data.csv") %>%
  mutate(Allele = factor(Allele,
        levels = rev(c("1_1","1_2","2_1","2_2","3_1","3_2","4_1","4_2","5_1","5_2",
                   "6_1","6_2","7_1","8_1","8_2","9_1","9_2","10_1","10_2",
                   "11_1","11_2","12_1","12_2","13_1","13_2","14_1","14_2","15_1","15_2",
                   "16_1","16_2","17_1","18_1","18_2")), ordered = T)) %>%
  mutate(Population = factor(Population))

pdf("SX_ridgeline.pdf", height = 11, width = 8)
ggplot(data, aes(x = Position, y = Allele, height = Point, scale = .007, 
            color = Population, fill = Population)) +
  geom_ridgeline(alpha = 1, linewidth = .001, show.legend = FALSE) + 
  labs(x = "G-rich Spacer Number", y = "Length (Amino Acids)", subtitle = "Allele") +
  scale_x_continuous(breaks = seq(0,180,10), labels = c("0", sapply(seq(20,180,20), function(x) c("", x))), expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_color_manual(values = c("white","white")) +
  scale_fill_manual(values = c("navyblue","seagreen")) +
  theme_classic() +
  theme(plot.subtitle = element_text(hjust = -.08, vjust = -1, size = 15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(vjust = 0),
        axis.ticks.y = element_blank())
  ggsave("SX_ridgeline.png", height = 11, width = 8)
dev.off()
