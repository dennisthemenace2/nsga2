###nsga 2 R implementation
##binary encoding



BITS = 16

##init population
initP=function(N,p, boundmin, boundmax){

    
  P = list()
  for(i in 1:N){
    np = list()
    #np$params = runif(p, min = min, max = max)
    
    params = paste(sample(0:1,BITS*p,replace=T),sep='',collapse='' ) 
    
    np$params = params
    
    np$rank= 0

    P[[i]] = np 
  }
  P
}
BinToDec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}

decode=function(params, boundmin,  boundmax){
  par = c()
  p = nchar(params) / BITS
  for( i in 0:(p-1) ){
    v = BinToDec( substring(params, i*BITS +1 , i*BITS + BITS ) )
    
    sc =   boundmin + (( boundmax- boundmin)/((2.0**BITS)-1.0)) * v
    par = c(par, sc)
  }
  
  par
}

###evaluate error funcion
evaluate = function(P,fn, boundmin, boundmax){
  for(i in 1:length(P)){
    cp = P[[i]]
    ###decode params
    dparams = decode(cp$params, boundmin,  boundmax)
    
    cp$err = fn(dparams)
    P[[i]] = cp
  } 
  P
}


dominates =function(p,q){
  any((p$err-q$err)<0)
}

fastNonDominatedSort = function(P){
  
  firstfront = c()
  
  for( x in 1:length(P)){
    p = P[[x]]
    p$domcnt = 0
    p$domset = c()
    
    for(y in 1:length(P)){
      q=P[[y]]
      if(dominates(p,q)){
        p$domset = c(p$domset,y) ##append idx of dominated solution
      }else if(dominates(q,p)){
        p$domcnt = p$domcnt+1 ##increase counter for solutions which dominate us
      }
      
    #  P[[y]]= q
    }
    
    if(p$domcnt==0){
      p$rank = 1 
      firstfront = c(firstfront, x) ### append index of first front solution
    }
    
    P[[x]] = p
  }
  resfront = list() ### prepare result list
  resfront[[1]] = firstfront;
  
  nextfront = firstfront;
  
  while(length(nextfront)>0 ){

    addfront = c()
    
    for(x in 1:length(nextfront)){
      p = P[[ nextfront[x] ]] ##
      ##nor for all q
      if(length(p$domset) > 0){
        for(y in 1:length(p$domset) ){
          q = P[[ p$domset[y] ]]
          q$domcnt =  q$domcnt-1
          if(q$domcnt == 0){ ##member of next front
            q$rank = length(resfront) +1
            addfront = c(addfront,p$domset[y] )
          }
          
          P[[p$domset[y]]] = q
        }
      }

    }## end for

    resfront[[length(resfront) +1]] = addfront ##pushback result might be NULL but noll is not added in R nice
    nextfront = addfront ##set next front

  }

  list('P'=P,'fronts' =resfront)
}


crossover = function(p,q,cros){
  if(runif(1,0,1) < cros){ ##this should be removed
    return (p)
  }
  params = c()
  
  for(i in 1:nchar(p) ){
    if(runif(1,0,1) < 0.5){
      params=c(params, substr(p, i, i) )
    }else{
      params=c(params, substr(q, i, i) )
    }
  }
  
  paste (params, sep = "", collapse = '')
}

point_mutation =function(p,mutation){
  params = c()
  for(i in 1:nchar(p) ){
    if(runif(1,0,1)<mutation){
      b = as.numeric(substr(p, i, i))
      if(b==1){
        params=c(params,0  )  
      }else{
        params=c(params,1  )
      }
    }else{
      params=c(params, substr(p, i, i) )
    }
    
  }
  
  paste (params, sep = "", collapse = '')
}

reproduce= function(P,selected, popsize, cros,mutation){
  childrens = list()
  
  while(length(childrens)<popsize){
    idx = sample(1:length(selected),2,replace=T )
    
    p = P[[ selected[idx[1]] ]]
    q = P[[ selected[idx[2]] ]]
    
    child = list()
    child$params = crossover( p$params, q$params , cros )
    child$params = point_mutation(child$params,mutation)
    
    
    childrens[[length(childrens)+1]] = child;
  }
  
 childrens
  
}

calculate_crowding_distance = function(union,cf,objectives){
  
  if(length(cf)<=2){
    union[[cf[1]]]$dist = Inf
    if(length(cf)==1){
      return (union)
    }
    
    union[[cf[2]]]$dist = Inf
    return (union)
  }
  
  for(i in 1:length(cf) ){ ##init distance measure to 0
    p = union[[cf[i]]]
    p$dist = 0
    union[[cf[i]]] = p
  }
  
  for(k in 1:objectives){
    
    values = c()
    for(i in 1:length(cf) ){
      p = union[[cf[i]]]
      values = c(values, p$err[k] )
    }
    sr = sort(values,index.return =T)
    union[[ cf[ sr$ix[1] ] ]]$dist = Inf
    union[[ cf[ sr$ix[length(sr$ix)] ] ]]$dist = Inf
    reg = max(values) - min(values)
    
    for(w in 2:(length(sr$ix)-1) ){
      union[[cf[sr$ix[w]]]]$dist= union[[cf[sr$ix[w]]]]$dist+ ( union[[cf[sr$ix[w+1]]]]$err[k] - union[[cf[sr$ix[w-1]]]]$err[k] )  /reg
    }
    
  }
 
  union
}


select_parents = function(union,fronts, popsize,objectives){
  
  parents =c()
    
  for(i in 1:length(fronts) ){
    cf = fronts[[i]]
    ##calc crowding distance
    union = calculate_crowding_distance(union,cf,objectives);
    
    if(length(parents)+ length(cf) > popsize){
      break;
    }
    parents = c(parents,cf)
  }
    
  ##pick up the solutions with highest front
  if(length(parents)< popsize){
    dif = popsize - length(parents)
    values = c()
    for(i in 1:length(cf)){
      values = c(values, union[[cf[i]]]$dist)
    }
    rs = sort(values,index.return =T, decreasing = TRUE)
    parents = c(parents, cf[ rs$ix[1:dif] ] )
  }
  
  ret = union[parents] # single brackets since these are more...
  ret
}

search = function(generations, popsize, cros,mutation,params,fn,objectives, boundmin=0, boundmax=1,tracDHV = FALSE  ){
 
  
  P = initP(popsize, params, boundmin, boundmax); # init population
  
  P = evaluate(P, fn, boundmin, boundmax)
  ##fast non dominated sort
  P = fastNonDominatedSort(P)$P
  
  a = sample(1:popsize,popsize ,replace=T )
  b = sample(1:popsize,popsize ,replace=T )
  selected =c()
  for( i in 1:popsize ){
    if( P[[a[i] ]]$rank < P[[b[i] ]]$rank ){
      selected = c(selected, a[i])
    }else{
      selected = c(selected, b[i])
    }
  }
  ##evaluate DHV
  if(tracDHV==T){
    require('mco') ##
    dhv = c()
  }
  
  children = reproduce(P,selected, popsize, cros,mutation)
  children= evaluate(children, fn, boundmin, boundmax)
  
  for( i in 1:generations){
 
    union = append( P , children)
    
    ret = fastNonDominatedSort(union)
    fronts = ret$fronts
    union = ret$P
    
    parents = select_parents(union,fronts, popsize,objectives)
    
    if(tracDHV==T){
      errmat = matrix( , nrow = 0, ncol = objectives )
      for(dv in 1:length(parents) ){
        #if(parents[[dv ]]$rank == 1){
          errmat = rbind(errmat, parents[[dv ]]$err)
      #  }
      }
      dhv = c(dhv,dominatedHypervolume(errmat))
    }
    
    
    
    a = sample(1:popsize,popsize ,replace=T )
    b = sample(1:popsize,popsize ,replace=T )
    selected =c()
    for( pop in 1:popsize ){
      if( parents[[a[pop] ]]$rank < parents[[b[pop] ]]$rank ){
        selected = c(selected, a[pop])
      }else if(parents[[a[pop] ]]$rank > parents[[b[pop] ]]$rank ){
        selected = c(selected, b[pop])
      }else{
        if(parents[[a[pop] ]]$dist > parents[[b[pop] ]]$dist){
          selected = c(selected, a[pop])
        }else{
          selected = c(selected, b[pop])
        }
      }
    }
    

    
    P=children
    children = reproduce(parents,selected, popsize, cros,mutation)
    children= evaluate(children, fn, boundmin, boundmax)
    
  }
  union = append( P , children)
  
  ret = fastNonDominatedSort(union)
  fronts = ret$fronts
  union = ret$P
  
  parents = select_parents(union,fronts, popsize,objectives)
  
  if(tracDHV==T){
    errmat = matrix( , nrow = 0, ncol = objectives )
    for(dv in 1:length(parents) ){
     #  if(parents[[dv ]]$rank == 1){
        errmat = rbind(errmat, parents[[dv ]]$err)
    #   }
    }
    dhv = c(dhv,dominatedHypervolume(errmat))
  }
  
  list(parents,dhv )
  
}

prob= function(params){
  err1 = sum(params^2)
  err2 = sum( (params-2)^2 )
  
  c(err1,err2)
}


res = search(100, 40,0.5,0.1,2,prob, 2,-10,10,T)
err = plotPopulation(res[[1]])

###plot population: // two objectives only...
plotPopulation=function(P){
  points =c()
  for(i in 1:length(P)){
    points = c(points, P[[i]]$err )
    
  }
  
 obj = length(P[[1]]$err)
 points  = matrix(points,nrow=obj)
 plot(points[1,],points[2,] , type='p')
 points
}

plotPopulation(P)

plotParameters=function(P,boundmin =-10,boundmax=10 ){
  points =c()
  for(i in 1:length(P)){
    points = c(points, decode( P[[i]]$params,boundmin,boundmax ) )
  }
  
  obj = 2 ##hard coded
  
  points  = matrix(points,nrow=obj)
  
  plot(points[1,],points[2,] , type='p')
  points
}

parents=res[[1]]
params = plotParameters(parents)
