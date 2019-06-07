require('Rcpp')

#Rcpp::sourceCpp('nsga2.cpp')
Rcpp::sourceCpp('nsga2RealEncoding.cpp')

###create a problem and call dat function!
GASettings <- list()
GASettings$parameterNr = 2
GASettings$objectiveNr = 2

GASettings$generationNr =100# GENRATIONS#12
GASettings$popSize = 40 #40

GASettings$cProb = 0.5  #crossover prob
GASettings$cDist = 5   #crossover distribution index
GASettings$mProb = 0.1  #mutation prob
GASettings$mDist = 10   #mutation distribution index

prob= function(params){
  #print(params)
  err1 = sum(params^2)
  err2 = sum( (params-2)^2 )
# mymoodModel$setParameters(params)
# sim = mymoodModel$simulate(LENGTH)
# err1 = sum ( (y$moodLevel - sim$moodLevel)^2 )
# err2 = sum ( (y$thoughts - sim$thoughts)^2 )
 
  c(err1,err2)
}

probVec = function(param){
  #print(param)
  t(apply(param, 1, prob)) 
}

# Start the clock!
ptm <- proc.time()

#result =  my_nsga2(prob,GASettings$parameterNr,GASettings$objectiveNr  ,
#                   GASettings$generationNr ,GASettings$popSize , GASettings$cProb,GASettings$cDist,
#                   GASettings$mProb,GASettings$mDist,-10, 10)

result =  my_nsga2(probVec,GASettings$parameterNr,GASettings$objectiveNr  ,
                   GASettings$generationNr ,GASettings$popSize , GASettings$cProb,GASettings$cDist,
                   GASettings$mProb,GASettings$mDist,0, 1)

###326.817


# Stop the clock
proc.time() - ptm

plot(result$error[,1],result$error[,2])

plot(result$parameters[,1],result$parameters[,2])


ptm <- proc.time()
r1 <- nsga2(probVec, 
            GASettings$parameterNr, 
            GASettings$objectiveNr,
            generations=GASettings$generationNr,
            popsize=GASettings$popSize, 
            cprob=GASettings$cProb,  
            cdist=GASettings$cDist,   
            mprob=GASettings$mProb,  
            mdist=GASettings$mDist,   
            lower.bounds=rep(0, GASettings$parameterNr),
            upper.bounds=rep(1, GASettings$parameterNr),
            vectorized=T)
# Stop the clock
proc.time() - ptm



plotStuff = function(res){
  
  v <- res$error[res$rank==1,]
  o <- res$rank==1
  d <- ncol(v)
  col <- ifelse(o, "red", "blue")
  pch <- ifelse(o, 4, 19)
  
  plot(v, col=col, pch=pch)
  
  v <- res$error
  
  ov <- v[o,]
  ov <- ov[order(ov[,1]),]
  lines (ov, col="red", type="s")
  
}


##multi threaded evaluation function

resultMT =  my_nsga2MT(prob,GASettings$parameterNr,GASettings$objectiveNr  ,
                              GASettings$generationNr ,GASettings$popSize , GASettings$cProb,
                              GASettings$mProb,-10, 10,1)

probVec = function(param){
#print(param)
  t(apply(param, 1, prob)) 
}

probVecPar = function(param){
  splitList <- split(param, 1:NROW(param))
  matrix(unlist(mclapply( splitList, prob )),ncol=ncol(param),byrow = TRUE ) ###ncol param is bug only owrd when  p = o
}

library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)


probVecDoPar = function(param){
  
 res <- foreach(i = 1:nrow(param),.combine=rbind,.export='prob') %dopar% {    
          prob(param[i,])
        }
 t(res)
}

#stopCluster(cl)
dhv =c()
reportFunction = function(error){
#  print('reported error')
#  print(error)
 # print(dominatedHypervolume(error))
  dhv<<- c(dhv,dominatedHypervolume(error))
}

ptm <- proc.time()
result = my_nsga2Vec(probVec,GASettings$parameterNr,GASettings$objectiveNr  ,
            GASettings$generationNr ,GASettings$popSize , GASettings$cProb,
            GASettings$mProb,-10, 10,reportFunction)
proc.time() - ptm

plot(1:length(dhv),dhv,type='o')


plot(result$error[,1],result$error[,2])




