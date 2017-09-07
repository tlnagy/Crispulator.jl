#setwd("~/Sandbox/Crispulator.jl/sim_data/")



count.to.matrix <- function( dir.name ) {
  N <- as.integer(strsplit(dir.name, "_")[[1]][2])
  fmt <- file.path(dir.name,"bc_count_%s.csv")
  
  df.count <- list()
  
  for(idx in 1:N) {
    df.count[[idx]] <- read.csv(sprintf(fmt, idx))
  }
  
  for(i in 1:(N-1)) {
    for(j in (i+1):N) {
      print(c(i, j, sum(df.count[[i]]$class != df.count[[j]]$class)))
    }
  }
  
  tmp <- df.count[[1]]
  
  df.label <- data.frame(gene=sprintf("G%04d",tmp$gene), sgRNA=sprintf("G%04d_%02d", tmp$gene, (tmp$barcodeid-1)%%10+1), class=tmp$class)
  df.low <- list()
  df.high <- list()
  
  for(i in 1:N) {
    tmp <- df.count[[i]]
    low.id <- sprintf("L%d", i)
    high.id <- sprintf("H%d", i)
    df.low[[low.id]] <- as.integer(tmp$counts_bin1-0.5)
    df.high[[high.id]] <- as.integer(tmp$counts_bin2-0.5)
  }
  
  df.final <- cbind(df.label, df.low, df.high)
  df.final
}

setwd("~/Sandbox/Crispulator.jl/simulation")
scenarios <- Sys.glob("data/*")

library(stringr)
for(dir.name in scenarios) {
  df.out <- count.to.matrix(dir.name)
  out.name <- paste0(str_replace(dir.name, "data", "matrix"), ".csv")
  write.csv(file=out.name, df.out)
}
