#Melanie Weilert
#January 2021
#Purpose: Define variables and filepaths relevant for analysis

tasks <- c('I_WT_D6CM','I_WT_S3MN')
color.vec<-c('#476c89','#914236')
names(color.vec)<-tasks

#Contribution filepaths
contrib.profile.bws<-lapply(tasks, function(x){
  paste0('preds/seq_width1000-lr0.001-lambda100-n_dil_layers9-conv_kernel_size7-tconv_kernel_size7-filters64/bw/',
         x, '.contrib.profile.bw')
})
names(contrib.profile.bws)<-tasks

contrib.counts.bws<-lapply(tasks, function(x){
  paste0('preds/seq_width1000-lr0.001-lambda100-n_dil_layers9-conv_kernel_size7-tconv_kernel_size7-filters64/bw/',
         x, '.contrib.counts.bw')
})
names(contrib.counts.bws)<-tasks

# Actual ChIP-nexus filepaths
actual.norm.bws<-list(
  I_WT_D6CM = 'data/bw/combined/I_wt_d6cm_log2_norm.bw',
  I_WT_S3MN = 'data/bw/combined/I_wt_s3mn_log2_norm.bw')

actual.bws<-list(
  I_WT_D6CM = 'data/bw/combined/I_wt_d6cm.bw',
  I_WT_S3MN = 'data/bw/combined/I_wt_s3mn.bw')

# Predicted ChIP-nexus filepaths
pred.bws<-lapply(tasks, function(x){
  list(paste0('preds/seq_width1000-lr0.001-lambda100-n_dil_layers9-conv_kernel_size7-tconv_kernel_size7-filters64/bw/_',
                 x,'.preds.pos.bw'))
})
names(pred.bws)<-tasks
