library(tidyverse)
library(RColorBrewer)



if(file.exists("benchmarks.txt")){
  data_info = read_tsv("benchmarks.txt")
}else{
  data_info = read_tsv("/home/projects/cu_10047/data/analyzed_samples/project_hierarchy/research_projects/CNV_benchmark/concordance/benchmarks.txt")
  
}



palette = c(brewer.pal(n=12, 'Paired')[1],brewer.pal(n=12, 'Paired')[5])

#Plot1. Tool vs memory
data_info %>% filter(!peak_memory =="." ) %>% ggplot() +
 geom_bar(aes(x=tool,y=as.numeric(peak_memory),fill=library),col="black",stat = "identity", position = "dodge") +
  xlab("Tool") + ylab("Peak Memory (Mb)") + theme_bw() + scale_fill_manual(values = palette) +
   theme(axis.text.x = element_text(angle = 90)) -> plot1


#Plot2. Tool vs Run time
data_info %>% filter(!cputime =="." ) %>% ggplot() +
  geom_bar(aes(x=tool,y=as.numeric(cputime)/60,fill=library),col="black",stat = "identity", position = "dodge") +
  xlab("Tool") + ylab("CPU time (h) in a 28-core machine") + theme_bw() + scale_fill_manual(values = palette) +
   theme(axis.text.x = element_text(angle = 90)) -> plot2



require(gridExtra)
png("Benchmark_time_memory.png") 
grid.arrange(plot1 + theme(legend.position = "none"), plot2, ncol=2)
dev.off() 
pdf("Benchmark_time_memory.pdf") 
grid.arrange(plot1 + theme(legend.position = "none"), plot2, ncol=2)
dev.off() 