source("CRCC_Analysis/040_Extract_results_Import_data.R")
source("CRCC_Analysis/041_Extract_results_functions.R")

out <- readRDS("CRCC_Analysis/Output/RDS/All_runs_Poseidon/Best_Poseidon.RDS")

# extract densities -------------------------------------------------------

clcl <- POSEIDON::estimate_clust_poseidon(out)

for(kk in 1:length(unique(clcl$col_cl))){
  cat(kk)
  G <- gg_draw_dens_poseidon(res = out,pixels = pixels,k = kk)
  #G+ggview::canvas(h=5,w=20)
  ggsave(plot = G, filename = paste0("CRCC_Analysis/Output/PLOT/Gof_density/CC-",kk,".eps"),
         height = 5,width = 20,device = cairo_ps)
  ggsave(plot = G, filename = paste0("CRCC_Analysis/Output/PLOT/Gof_density/CC-",kk,".pdf"),height = 5,width = 20)
  ggsave(plot = G, filename = paste0("CRCC_Analysis/Output/PLOT/Gof_density/CC-",kk,".png"),height = 5,width = 20)
}




# extract dots -------------------------------------------------------

clcl <- POSEIDON::estimate_clust_poseidon(out)

allx <- list(allx__g,allx__l,allx__p)
limx <- list(limx__g,limx__l,limx__p)

for(kk in 1:length(unique(clcl$col_cl))){
  cat(kk)

  G <- gg_draw_dots(out = out, pixels = pixels,ind = kk,median = TRUE,
                    allx =  allx, limx =  limx)
  #G+ggview::canvas(h=5,w=20)
  ggsave(plot = G, filename = paste0("CRCC_Analysis/Output/PLOT/Gof_density/CC-",kk,"dots_median.eps"),
         height = 5,width = 20,device = cairo_ps)
  ggsave(plot = G, filename = paste0("CRCC_Analysis/Output/PLOT/Gof_density/CC-",kk,"dots_median.pdf"),height = 5,width = 20)
  ggsave(plot = G, filename = paste0("CRCC_Analysis/Output/PLOT/Gof_density/CC-",kk,"dots_median.png"),height = 5,width = 20)
}


for(kk in 1:length(unique(clcl$col_cl))){
  cat(kk)
  G <- gg_draw_dots(out = out, pixels = pixels,ind = kk,median = FALSE,
                    allx =  allx, limx =  limx)
  #G+ggview::canvas(h=5,w=20)
  ggsave(plot = G, filename = paste0("CRCC_Analysis/Output/PLOT/Gof_density/CC-",kk,"dots_mean.eps"),
         height = 5,width = 20,device = cairo_ps)
  ggsave(plot = G, filename = paste0("CRCC_Analysis/Output/PLOT/Gof_density/CC-",kk,"dots_mean.pdf"),height = 5,width = 20)
  ggsave(plot = G, filename = paste0("CRCC_Analysis/Output/PLOT/Gof_density/CC-",kk,"dots_mean.png"),height = 5,width = 20)
}



for(kk in 1:length(unique(clcl$col_cl))){
  cat(kk)
  G <- gg_draw_dots_detailed_axis(out = out, pixels = pixels,ind = kk,median = TRUE,
                    allx =  allx, limx =  limx)
  #G+ggview::canvas(h=5,w=20)
  ggsave(plot = G, filename = paste0("CRCC_Analysis/Output/PLOT/Gof_density/CC-",kk,"dots_median_detailed.eps"),
         height = 5,width = 30,device = cairo_ps)
  ggsave(plot = G, filename = paste0("CRCC_Analysis/Output/PLOT/Gof_density/CC-",kk,"dots_median_detailed.pdf"),height = 5,width = 30)
  ggsave(plot = G, filename = paste0("CRCC_Analysis/Output/PLOT/Gof_density/CC-",kk,"dots_median_detailed.png"),height = 5,width = 30)
}
