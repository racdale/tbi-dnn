se = function(x) 1.96 * sd(x) / sqrt(length(x))
targetLayer = 10

profile = function(targetLayer, self=c(T,F)) {
  s_mean = aggregate(residH~rel_dist+tb, data=processedData[processedData$layer%in%targetLayer&processedData$self%in%self,], mean)
  print(s_mean)
  s_se = aggregate(residH~rel_dist+tb, data=processedData[processedData$layer%in%targetLayer&processedData$self%in%self,], se)
  dat = merge(s_mean, s_se, by=c("rel_dist", "tb"))
  colnames(dat) = c("rel_dist", "tb", "residH_mean", "residH_se")
  
  # Residualized H x relative distance (k)
  p1 = ggplot(dat, aes(x = rel_dist, y = residH_mean)) +
    geom_ribbon(aes(ymin = residH_mean - residH_se, ymax = residH_mean + residH_se, fill = as.factor(tb)), alpha = 0.3) +
    geom_line(aes(color = as.factor(tb))) +
    labs(x = "Relative distance (k)", y = "Residualized entropy", color = "Group", fill = "Group") +
    scale_color_manual(values = c("red", "blue"), labels = c("Control", "TBI")) +
    scale_fill_manual(values = c("red", "blue"), labels = c("Control", "TBI")) +
    labs(title=paste('Layer',targetLayer))+theme(legend.position = "none")
  if (targetLayer[1]==1) {
    p1=p1+geom_label(aes(label = "Control"), 
                     x = -8, y = -.082, hjust = 0, vjust = 0, color='red')
    p1=p1+geom_label(aes(label = "TBI"), 
                     x = -8, y = -.085, hjust = 0, vjust = 0, color='blue')
  }
  return(p1)
}

