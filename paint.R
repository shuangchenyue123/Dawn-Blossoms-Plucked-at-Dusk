library(ggplot2)
data = read.csv("Events.tsv",sep = "\t")
range1 = max(data$gene_present_time.mya.)
range2 = max(data$horizontal_transfer_time.mya.)
jitter = (range1 + range2)/2/100
p1 <- ggplot(data=data, aes(x=gene_present_time.mya., y=horizontal_transfer_time.mya.,alpha=0.1)) +
  geom_jitter(aes(size=event_confidence,color=cluster), width = jitter, height = jitter)+
  scale_size_continuous(range = c(1,3)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black")) +
  guides(size=guide_legend(order=1))+
  guides(alpha = "none")+
  xlab("gene_present_time(mya)")+
  ylab("horizontal_transfer_time(mya)")


p2 <- ggplot(data=data, aes(x=gene_present_time.mya., y=horizontal_transfer_time.mya.,alpha=0.1)) +
  geom_jitter(aes(size=event_confidence,color=orthogroup), width = jitter, height = jitter)+
  scale_size_continuous(range = c(1,3)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black")) +
  guides(size=guide_legend(order=1))+
  guides(alpha = "none", color = "none")+
  xlab("gene_present_time(mya)")+
  ylab("horizontal_transfer_time(mya)")

ggsave("Event_cluster.png", p1, width = 10, height = 8, dpi=1200)
ggsave("Event_ortho.png", p2, width = 10, height = 8, dpi=1200)