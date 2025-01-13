library (ggplot2)
setwd("C:/Users/Aurelio/Desktop/Organoidi_2024/")

timepoints <- c("23days","1month","1.5month","2month","3month","4month","5month","6month")

name_new_dir <- paste(getwd(),"/Results/","UMAP_dal_paper",sep="")
if (dir.exists(name_new_dir)==FALSE) {
  dir.create(name_new_dir)
} 

recreate_UMAP <- function(name) {
  UMAP_coordinates <- read.delim(paste(getwd(),"/Data/cluster/umap_",name,".txt",sep=""))
  UMAP_coordinates <- UMAP_coordinates[-1,]
  for (i in 2:7) {if (i!=4) {UMAP_coordinates[,i]<-as.numeric(UMAP_coordinates[,i])}}
  
  UMAP_plot <- ggplot(UMAP_coordinates, aes(x = X, y = Y, color = CellType)) +
    geom_point(size = 1, alpha = 0.8) +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 24),
      legend.text = element_text(size = 20)
    ) +
    guides(color = guide_legend(override.aes = list(size = 5))) # Cambia la dimensione dei punti nella legenda

  jpeg(paste(name_new_dir,"/",name,".jpg",sep=""), width = 1920, height = 1080)
  print(UMAP_plot)
  dev.off()
}

for (i in 1:length(timepoints)) {
  recreate_UMAP(timepoints[i])
  print(timepoints[i])
}

