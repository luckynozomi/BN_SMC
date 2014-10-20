args<- commandArgs(trailingOnly = TRUE)

TOTAL_NO_DATASET<-as.numeric(args[1])
NO_NODES<-as.numeric(args[2])
NO_OBS<-as.numeric(args[3])
MAX_ITERATION<-as.numeric(args[4])
NO_CHAINS<-as.numeric(args[5])
CUTOFF<-as.numeric(args[6])
NO_CORES<-as.numeric(args[7])
TEMPERATURE<-as.numeric(args[8])
DATA_NAME<-args[9]

for(i in 1:TOTAL_NO_DATASET)
{
  DATA_PATH<-sprintf("./data/%s/%s/%s_data_%d.txt",args[3],args[9],args[9],i)
  PARAM_PATH<-sprintf("./data/%s/%s/%s_param_%d.txt",args[3],args[9],args[9],i)
  OUTPUT_PATH<-sprintf("./output/%s/%s_%d_ch%d",args[3],args[9],i,NO_CHAINS)
  system(paste(sprintf("./a %d %s %d %d %f %d %f %s %s %s",NO_NODES,args[3],
            MAX_ITERATION,NO_CHAINS,CUTOFF,NO_CORES,TEMPERATURE,DATA_PATH,PARAM_PATH,OUTPUT_PATH)),
         ignore.stdout = T,ignore.stderr = T)  
  
}





