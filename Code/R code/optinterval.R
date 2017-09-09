library("cairoDevice")

#####constant fitter

constantobjective = function(N,A){
  value = N*(log(N/A))-N
  value[is.nan(value)]<-0
  return(value)
}

##############################exponential fitter
powfunc= function(x,t){
  a=x[1]
  b=x[2]
  c=x[3]
  tn=t[length(t)]
  eb=exp(1)^b
  sum(log(a*t^eb+c))-((a*tn^(eb+1))/(eb+1)+c*tn)
}

powgunc= function(x,t){
  a=x[1]
  b=x[2]
  c=x[3]
  eb=exp(1)^b
  atbc=a*t^(eb)+c
  teb=t^eb
  
  tn=t[length(t)]
  logt=log(t)
  fa=sum(teb/atbc)-((tn^(eb+1))/(eb+1))
  fb=sum((a*eb*t^(eb+1)*(eb*logt+logt-1))/(eb+1)^2)-(a*eb*tn^(eb+1)*((eb+1)*log(tn)-1))/(eb+1)^2
  fc=sum(1/atbc)-tn
  return(c(fa,fb,fc))
}


powfit <- function(N,x){
  if(length(x)<15){
    return(list("a"=0,"b"=0,"c"=length(x)/x[length(x)],"direction"="forward","cost"=constantobjective(length(N),x[length(x)])))
  }
  x=rep(x,N)
  
  #fit the data normally
  fit <- optim(par=c(1,1,length(x)/x[length(x)]),fn=powfunc,gr=powgunc,control = list(fnscale = -1),t=x)
  pars=fit$par
  a=pars[[1]]
  b=exp(1)^pars[[2]]
  c=pars[[3]]
  cost=fit$value
  forward=list("a"=a,"b"=b,"c"=c,"direction"="forward", "cost"=cost)
  if(sum((a*c(0,x)^b+c)<0)>0){
   forward$cost=-Inf
  }
  
  #fit the data in reverse
  x=x[length(x)]-rev(x[1:length(x)-1])
  fit <- optim(par=c(1,1,length(x)/x[length(x)]),fn=powfunc,gr=powgunc,control = list(fnscale = -1),t=x)
  pars=fit$par
  a=pars[[1]]
  b=exp(1)^pars[[2]]
  c=pars[[3]]
  cost=fit$value
  reverse=list("a"=a,"b"=b,"c"=c,"direction" = "reverse","cost"=cost)
  if(sum((a*c(0,x)^b+c)<0)>0){
    reverse$cost=-Inf
  }
  
  if (forward$cost>reverse$cost){
    return(forward)
  }else{
    return(reverse)
  }
}
###########################FRED fitter

FREDfunc= function(x,t){
  A=exp(1)^x[1]
  lambda=exp(1)^x[2]
  tau1=exp(1)^x[3]
  tau2=exp(1)^x[4]
  b=exp(1)^x[5]
  sum(log(A*lmabda*exp(1)^(-tau1/t-t/tau2)))-(2*a*lambda*(4*tau1*tau2)^0.5)*besselI((4*tau1/tau2)^0.5,1)-b*t[length(t)]
}

FREDfit <- function(N,x){
  if(length(x)<5){
    return(list("A"=0,"lambda"=0,"tau1"=0,"tau2"=0,"b"=0,"cost"=-Inf))
  }
  x=rep(x,N)
  
  #fit the data normally
  fit <- optim(par=c(1,1,1,1,log(length(x)/x[length(x)])),fn=powfunc,control = list(fnscale = -1),t=x)
  pars=fit$par
  A=exp(1)^pars[[1]]
  lambda=exp(1)^pars[[2]]
  tau1=exp(1)^pars[[3]]
  tau2=exp(1)^pars[[4]]
  b=exp(1)^pars[[5]]
  cost=fit$value
  return(list("A"=A,"lambda"=lambda,"tau1"=tau1,"tau2"=tau2,"b"=b,"cost"=cost))
}

##############################exponential fitter
exp.f.iwls= function(x,Yk,xk){
  a=x[1]
  b=x[2]
  c=x[3]
  f=a*exp(1)^(b*xk)+c
  total=sum(1/f)
  sum((f-Yk)^2/(f*total))
}

exp.g.iwls= function(x,Yk,xk){
  a=x[1]
  b=x[2]
  c=x[3]
  fa=exp(1)^(b*xk)
  fb=a*xk*exp(1)^(b*xk)
  fc=1
  f=a*exp(1)^(b*xk)+c
  f2=f^2
  f3=f2*f
  y2=Yk^2
  total=sum(1/f)
  part1=total*(f2-y2)
  part2=(f3-2*y*f2+y2*f)
  denom=f2*total^2
  da=fa*part1-sum(fa/f2)*part2
  db=fb*part1-sum(fb/f2)*part2
  dc=fc*part1-sum(fc/f2)*part2
  return(c(da,db,dc)/denom)
}

exp.f.ML= function(x,t){
  a=x[1]
  b=x[2]
  c=x[3]
  tn=t[length(t)]
  eb=exp(1)^b
  ec=exp(1)^c
  sum(log(a*exp(1)^(eb*t)+ec))-(a/(eb)*(exp(1)^(eb*tn)-1)+ec*tn)
}

exp.g.ML= function(x,t){
  a=x[1]
  b=x[2]
  c=x[3]
  eb=exp(1)^b
  ec=exp(1)^c
  teb=exp(1)^(t*b)
  tn=t[length(t)]
  eteb=exp(1)^(t*eb)
  denom=(a*eteb+ec)
  etneb=exp(1)^(tn*eb)
  fa=sum(eteb/denom)-(1/eb)*(exp(1)^(eb*tn)-1)
  fb=sum((a*t*eb*eteb)/denom)-(a/eb)*(tn*eb*etneb-etneb+1)
  fc=sum(1/denom)-tn*ec
  return(c(fa,fb,fc))
}

expfit <- function(N,x){
  x=rep(x,N)
  if(length(x)<15){
    return(list("a"=0,"b"=0,"c"=length(x)/x[length(x)],"direction"="forward","cost"=constantobjective(length(N),x[length(x)])))
  }
  #fit the data normally
  fit <- optim(par=c(0,0,log(length(x)/x[length(x)])),fn=exp.f.ML,gr=exp.g.ML,control = list(fnscale = -1),t=x)
  pars=fit$par
  a=pars[[1]]
  b=exp(1)^pars[[2]]
  c=exp(1)^pars[[3]]
  cost=fit$value
  forward=list("a"=a,"b"=b,"c"=c,"direction"="forward", "cost"=cost)
  
  #fit the data in reverse
  x=x[length(x)]-rev(x[1:length(x)-1])
  fit <- optim(par=c(0,0,log(length(x)/x[length(x)])),fn=exp.f.ML,gr=exp.g.ML,control = list(fnscale = -1),t=x)
  pars=fit$par
  a=pars[[1]]
  b=exp(1)^pars[[2]]
  c=exp(1)^pars[[3]]
  #cost=sum(log(a*exp(1)^(b*x)+c))-((a/b)*(exp(1)^(b*x[length(x)])-1)+c*x[length(x)])
  cost=fit$value
  reverse=list("a"=a,"b"=b,"c"=c,"direction" = "reverse","cost"=cost)
  
  if (forward$cost>reverse$cost){
    return(forward)
  }else{
    return(reverse)
  }
  
}


expfit2 <- function(N,x){
  
  if(length(x)<20){
    return(list("a"=0,"b"=0,"c"=0,"direction"="forward","cost"=-Inf))
  }
  x=rep(x,N)
  bins=as.integer(max(2,min(length(x)/10,20)))
  
  # number of breaks needs to depend on the time length to give accurate estimates
  histogram = hist(x, breaks =seq(0,x[length(x)],length.out=bins+1), plot = FALSE)
  L = x[length(x)]         # total interval length
  N = length(histogram$counts)      # number of blocks on this interval
  binwidth=L/N
  k = 1:N                           # k'th block index (1<=k<=N)
  xk = (k - 1/2)*binwidth            # midpoints of blocks (from page 367)
  Yk = histogram$counts/binwidth   
  fit <- optim(par=c(1,0,10),fn=func,gr=gunc,Yk=Yk,xk=xk)
  pars=fit$par
  a=pars[[1]]
  b=pars[[2]]
  c=pars[[3]]
  if(sum((a*exp(1)^(b*x)+c)<0)>0){
    return(list("a"=0,"b"=0,"c"=0,"cost"=-Inf))
  }
  cost=sum(log(a*exp(1)^(b*x)+c))-((a/b)*exp(1)^(b*x[length(x)])+c*x[length(x)]-a/b)
  return(list("a"=a,"b"=b,"c"=c,"cost"=cost))
}

#### iwls linear fitter ###

iwls = function(N,t) {
  if (length(t)<40){ # if no data is recieved, set cost to 0
    return(list("a"=0,"b"=0,"cost"=-Inf))
  }  
  #if (length(t)<){ #if only one point is supplied, treat as constant intensity
  #  return(list("a"=0,"b"=N/t[1],"cost"=log(N/t[1])))
  #}
  x=rep(t,N)
  bins=as.integer(max(length(x)/7,2))
  # number of breaks needs to depend on the time length to give accurate estimates
  histogram = hist(x, breaks =seq(0,x[length(x)],length.out=bins+1), plot = FALSE)
  L = max(histogram$breaks)         # total interval length
  N = length(histogram$counts)      # number of blocks on this interval
  binwidth=L/N
  k = 1:N                           # k'th block index (1<=k<=N)
  xk = (k - 1/2)*binwidth            # midpoints of blocks (from page 367)
  Yk = histogram$counts/binwidth            # event counts in bins
  epsilon=1
  est.weights=rep(1,N)
  M=x[length(x)]
  S=M^2/2
  old.cost=-Inf
  while(epsilon>1e-4){
    m = lm(Yk ~ xk,weights=est.weights)
    est.weights = 1/m$fitted.values
    coef=m$coefficients
    new.cost=(sum(log(coef[[2]]*x+coef[[1]]))-coef[[2]]*S-coef[[1]]*M)
    if (is.nan(new.cost) || sum(est.weights<0)>0){
      return(list("a"=0,"b"=0,"cost"=-Inf))
    }
    epsilon = new.cost-old.cost
    old.cost=new.cost
  }
  
  list("a"=coef[2],"b"=coef[1],"cost"=new.cost)
}

#######linear fitter

linfit <-function(N,t){
  
  if (length(t)==0){ # if no data is recieved, set cost to 0
    return(list("a"=0,"b"=0,"cost"=-Inf))
  }  
  if (length(t)==1){ #if only one point is supplied, treat as constant intensity
    return(list("a"=0,"b"=N/t[1],"cost"=log(N/t[1])))
  }
  
  t=rep(t,N)
  
  epsilon  <- 1
  n        <- length(t)
  i        <- 1
  M        <- t[n]
  S        <- t[n]^2/2
  
  #start with some initial values for a and b
  coef     <- matrix(c(0,length(t)/M),nrow=2)
  
  new.cost <- -Inf
  
  while (abs(epsilon) >10e-4 && i<100){
    old.cost <- new.cost
    
    #first derivatives
    f=matrix(c(sum(t/(coef[1]*t+coef[2]))-S,sum(1/(coef[1]*t+coef[2]))-M),nrow=2)
    
    #hessian values
    fa <- -sum((t/(coef[1]*t+coef[2]))^2)
    fb <- -sum(t/(coef[1]*t+coef[2])^2)
    ga <- fb
    gb <- -sum(1/(coef[1]*t+coef[2])^2)
    
    #create the inverse of the hessian
    invhess <- 1/(fa*gb-fb*ga)*matrix(c(gb,-fb,-ga,fa),nrow=2,ncol=2,byrow=TRUE) 
    
    #run newtons method
    coef=coef-invhess%*%f 
    
    if(coef[2]<0 || (M*coef[1]+coef[2])<0){#if best line has a negative intensity, dont use linear
      return(list("a"=0,"b"=0,"cost"=-Inf))
    }
    #calculate the new cost
    new.cost <- sum(log(coef[1]*t+coef[2]))-coef[1]*S-coef[2]*M
    epsilon  <- new.cost-old.cost
    
    i <- i+1
  }
  list("a"=coef[1],"b"=coef[2],"cost"=new.cost)
}

#S3 generic plot function 

"plot.BB" <- function(x,show=c("hist","blocks"),binwidth=NULL,bins=NULL,ylim=NULL,xlim=NULL,xact=NULL,main="Bayesian Blocks",xlab="Time",ylab="Intensity") {
  
  data       <- x$data
  n          <- length(data)
  intensity  <- round(x$N/x$A,4)
  legend     <- vector()
  legend.col <- vector()
  
  #check if xlim is supplied
  if (!is.null(xlim)){
    lowerbound=xlim[1]
    upperbound=xlim[2]
  }
  else{
    lowerbound <- data[1]-(data[2]-data[1])/2
    upperbound <- data[n]+(data[n]-data[n-1])/2  
  }
  
  #check if ylim is supplied
  if(is.null(ylim)){
    ylim=c(0,max(intensity)*0.75)
  }
  
  #check if the hist bins are supplied
  if (is.null(binwidth) & is.null(bins)){
    binwidth=data[n]/100
    bins=100
  }
  else if (is.null(bins)){
    bins   <- round(data[n]/binwidth,0)
  }
  else{
    binwidth=data[n]/bins
  }
  
  
  #initialize plot
  plot(NA, ylim = ylim, xlim = c(lowerbound,upperbound),xaxt=xact,bty="n",
       main=main,xlab=xlab,ylab=ylab)  
  grid(NA,NULL)#add grid for y-axis
  
  
  
  #plot histogram if requested
  if ("hist" %in% show){
    histdata <- rep(data,x$N)
    end      <- histdata[length(histdata)]
    i        <- 0:(end/binwidth)
    
    height   <- hist(histdata, breaks=c(i,bins+1)*binwidth,plot=FALSE)$counts/binwidth
    rect(xleft=i*binwidth,ybottom=0,xright= (i+1)*binwidth, ytop=height,col="grey70",border="grey40")
    legend=c(legend,"Binned Data")
    legend.col=c(legend.col,"grey70")
  }
  
  
  #plot individual intensities if requested
  if ("points" %in% show){ 
    segments(x0=c(0,data[1:(length(data)-1)]), y0=intensity[1:length(intensity)], 
             x1=data[1:length(data)], y1=intensity[1:length(intensity)],
             col = "cornflowerblue", lwd = 3)
    legend=c(legend,"Data Points")
    legend.col=c(legend.col,"cornflowerblue")
  }
  
  
  #plot blocks if requested
  if ("blocks" %in% show){
    #plot constant blocks
    indices <- which(x$type=="constant")
    if (length(indices)>0){
      b=sapply( x$params [indices], "[[" , "b" )
      segments(x0=data[x$left[indices]], y0=b, 
               x1=data[x$right[indices]], y1=b,
               col = "red", lwd = 3)
    }
    #plot linear blocks
    indices <- which(x$type=="linear")
    if(length(indices)>0){
      a=sapply( x$params[indices], "[[" , "a" )
      b=sapply( x$params[indices], "[[" , "b" )
      lefts=data[x$left[indices]]
      rights=data[x$right[indices]]
      segments(x0=lefts, y0=b, 
               x1=rights, y1=a*(rights-lefts)+b,
               col = "red", lwd = 3)
    }
    
    indices <- which(x$type=="exponential")
    if(length(indices)>0){
      a=sapply( x$params[indices], "[[" , "a" )
      b=sapply( x$params[indices], "[[" , "b" )
      c=sapply( x$params[indices], "[[" , "c" )
      dir=sapply( x$params[indices], "[[" , "direction" )
      lefts=data[x$left[indices]]
      rights=data[x$right[indices]]
      
      for (i in 1:length(a)){
        xv=seq(lefts[i],rights[i],length.out=100)
        if(dir[i]=="forward"){
          lines(xv,a[i]*exp(1)^(b[i]*(xv-lefts[i]))+c[i],col="red",lwd=3)
        }else{
          lines(rev(xv),a[i]*exp(1)^(b[i]*(rights[i]-rev(xv)))+c[i],col="red",lwd=3)
        }
      }
    }
    
    indices <- which(x$type=="power")
    if(length(indices)>0){
      a=sapply( x$params[indices], "[[" , "a" )
      b=sapply( x$params[indices], "[[" , "b" )
      c=sapply( x$params[indices], "[[" , "c" )
      dir=sapply( x$params[indices], "[[" , "direction" )
      lefts=data[x$left[indices]]
      rights=data[x$right[indices]]
      
      for (i in 1:length(a)){
        xv=seq(lefts[i],rights[i],length.out=100)
        if(dir[i]=="forward"){
          lines(xv,a[i]*(xv-lefts[i])^b[i]+c[i],col="red",lwd=3)
        }else{
          lines(rev(xv),a[i]*(rights[i]-rev(xv))^b[i]+c[i],col="red",lwd=3)
        }
      }
        
    }
    #plot FRED blocks
    indices <- which(x$type=="FRED")
    if(length(indices)>0){
      A=sapply( x$params[indices], "[[" , "A" )
      lambda=sapply( x$params[indices], "[[" , "lambda" )
      tau1=sapply( x$params[indices], "[[" , "tau1" )
      tau2=sapply( x$params[indices], "[[" , "tau2" )
      b=sapply( x$params[indices], "[[" , "b" )
      lefts=data[x$left[indices]]
      rights=data[x$right[indices]]
      for (i in 1:length(indices)){
        xv=seq(lefts[i],rights[i],length.out=100)
        tb=xv-lefts[i]
        lines(xv,a[i]*lambda[i]*exp(1)^(-tau1/tb-tb/tau2)+b[i],col="red",lwd=3)
      }
    }
    
    
    
    legend=c(legend,"Blocks")
    legend.col=c(legend.col,"red")
    
  }
  
  #plot legend
  legend(x=x$data[1],y=ylim[2]*0.98,legend= legend, lty=rep(1,length(legend)), 
         lwd=rep(3,length(legend)), col=legend.col) 
}


#S3 Generic summary function

summary.BB <- function(x) {
  cat(length(x$data), 'data points\n\n',
      length(x$left),'total blocks:\n\t',
      sum(x$type=='constant'),'constant blocks\n\t',
      sum(x$type=='linear'),'linear blocks\n\t',  
      sum(x$type=='exponential'),'exponential blocks\n\t',
      sum(x$type=='power'),'power function blocks\n\t',
      sum(x$type=='FRED'),'FRED blocks\n\n',
      floor((length(x$pruned)*100)/length(x$data)),'% of the points were pruned\n\n')
}


# OPTINTERVAL FUNCTION 
# An O(N^2) algorithm for finding the optimal partition of N data points on an interval
# INPUT:
# data is the vector ofcells
# N is the number of "observations" per point. This is the same "N" we've been using all along
# c is the block penalty term
# OUTPUT:
# a list of objects pertaining to the "BB" object

optinterval = function(data,N,c,type=c("constant"),pvalue=0.05,verbose=FALSE){
  start.time <- Sys.time()
  
  chi.lin    <- qchisq(1-pvalue,df=1)/2
  chi.pow    <- qchisq(1-pvalue,df=2)/2
  chi.exp    <- qchisq(1-pvalue,df=2)/2
  if("power" %in% type | "exponential" %in% type){
    maxpen   <- chi.pow +c  
  }else if("linear" %in% type){
    maxpen   <- chi.lin+c
  }else{
    maxpen   <- c
  }
  
  n          <- length(N)
  percent    <- floor(n/100)
  lowerbound <- data[1]-(data[2]-data[1])/2
  upperbound <- data[n]+(data[n]-data[n-1])/2
  xx         <- c(lowerbound,(data[1:(length(data)-1)] + data[2:(length(data))])/2,upperbound) # voronoi cell vertices
  A          <- diff(xx) # length of each voronoi cell
  
  opt           <- rep(0,n+1)
  lastchange    <- rep(1,n)
  changeA       <- rep(0,n) 
  changeN       <- rep(0,n)
  optint        <- matrix(-Inf,nrow=n,ncol=length(type))
  last          <- rep(0,length(type))
  lasttype      <- rep("None",n)
  lastparams    <- list()
  unpruned <- NULL
  endobj   <- rep(0,n)

  
  
  #begin looping through each point
  for (i in 1:n){
    
    unpruned          <- c(unpruned,i)
    
    ##### constant blocks ##########
    if ("constant" %in% type){
      changeA[unpruned]  <- changeA[unpruned] + A[i]
      changeN[unpruned]  <- changeN[unpruned] + N[i]
      optint[unpruned,which(type=="constant")]  <- opt[unpruned] + constantobjective(changeN[unpruned],changeA[unpruned])-c
      last[which(type=="constant")]             <- which.max(optint[unpruned,which(type=="constant")])
    }
    ################################
    
    ##### linear blocks ############
    if ("linear" %in% type){
      linblocks=list()
      for(j in unpruned){
        x=data[j:i]-max(data[j-1],0)
        linblocks[[j]]   <- linfit(N[j:i],x)
        #linblocks[[j]]   <- iwls(N[j:i],x)
        optint[j,which(type=="linear")]       <- opt[j] + linblocks[[j]][["cost"]]-c-chi.lin
      }
      last[which(type=="linear")]              <- which.max(optint[unpruned,which(type=="linear")])
    }
    
    ##### exp blocks ############
    if ("exponential" %in% type){
      expblocks=list()
      #print(unpruned)
      for(j in unpruned){
        x=data[j:i]-max(data[j-1],0)
        expblocks[[j]]   <- expfit(N[j:i],x)
        optint[j,which(type=="exponential")]       <- opt[j] + expblocks[[j]][["cost"]]-c-chi.exp
      }
      last[which(type=="exponential")]                  <- which.max(optint[unpruned,which(type=="exponential")])
    }
      ################################
    
    ##### power blocks ############
    if ("power" %in% type){
      powblocks=list()
      for(j in unpruned){
        x=data[j:i]-max(data[j-1],0)
        powblocks[[j]]   <- powfit(N[j:i],x)
        optint[j,which(type=="power")]       <- opt[j] + powblocks[[j]][["cost"]]-c-chi.pow
      }
      last[which(type=="power")]               <- which.max(optint[unpruned,which(type=="power")])
    }
    ################################  
    
    ##### FRED blocks ############
    if ("FRED" %in% type){
      FREDblocks=list()
      for(j in unpruned){
        x=data[j:i]-max(data[j-1],0)
        FREDblocks[[j]]   <- FREDfit(N[j:i],x)
        optint[j,which(type=="FRED")]       <- opt[j] + FREDblocks[[j]][["cost"]]-c-chi.exp
      }
      last[which(type=="FRED")]              <- which.max(optint[unpruned,which(type=="FRED")])
    }
    ################################
    
    ################################
    
    bestshape<-which.max(apply(optint,2,max))
    
    lastchange[i]     <- unpruned[bestshape]
    if(("FRED" %in% type) && (bestshape==which(type=="FRED"))){ #pow block is best
      lasttype[i]       <- "FRED"
      lastparams[[i]]   <- FREDblocks[[unpruned[which(type=="FRED")]]]
    }
    else if (("linear" %in% type) && (bestshape==which(type=="linear"))){#linear block is best
      lasttype[i]       <- "linear"
      lastparams[[i]]   <- linblocks[[unpruned[which(type=="linear")]]]
    }
    
    else if(("exponential" %in% type) && (bestshape==which(type=="exponential"))){ #exp block is best
      lasttype[i]       <- "exponential"
      lastparams[[i]]   <- expblocks[[unpruned[which(type=="exponential")]]]
    }
    
    else if(("power" %in% type) && (bestshape==which(type=="power"))){ #pow block is best
      lasttype[i]       <- "power"
      lastparams[[i]]   <- powblocks[[unpruned[which(type=="power")]]]
    }
    

    if(("constant" %in% type) && (bestshape==which(type=="constant"))){ #constant block is best
      lasttype[i]       <- "constant"
      lastparams[[i]]   <- list("b"= sum(N[lastchange[i]:i])/sum(A[lastchange[i]:i]))
    }
    
    opt[i+1]          <- max(optint[unpruned,bestshape])
    unpruned          <- unpruned[((optint[unpruned,bestshape]+maxpen-opt[i+1])>0)]

    if((verbose==TRUE) && (i %% percent==0)){#print out the progress of the algorithm
      cat("\n",round(100*i/n,0),"% of points completed",round(100*(i-length(unpruned))/i,0),"% pruned ")
    }
    
  }
  BBdata   <- list("data"         = data, 
                   "N"            = N,
                   "A"            = A,
                   "opt"          = opt,
                   "lastchange"   = lastchange,
                   "lastparams"   = lastparams,
                   "lasttype"     = lasttype,
                   "pruned"       = setdiff(1:i,unpruned))
  
  
  
  print(Sys.time()-start.time)
  model <- getStep(BBdata,n)
  summary(model)
  return(model)
}



#this function returns a BB object for the optimal blocks at any given cell location "n" where 1 < n < N
getStep <- function(BBobj,n){
  
  lasttype      <- BBobj$lasttype[1:n]
  lastchange    <- BBobj$lastchange[1:n]
  lastparams    <- BBobj$lastparams[1:n]
  
  lastchange<-lastchange-1
  
  left          <- vector()
  params        <- list()
  type          <- vector()
  type[1]       <- lasttype[n]
  left[1]       <- lastchange[n]
  params[[1]]   <- lastparams[[n]]
  i=1
  while (left[i] > 0) {
    left[i+1]      <- lastchange[left[i]]
    type[i+1]      <- lasttype[left[i]]
    params[[i+1]]  <- lastparams[[ left[i] ]]
    
    i <- i+1
  }
  left      <- rev(left)+1
  type      <- rev(type)
  params    <- rev(params)
  
  if (length(left) == 1){
    right <- n
  }
  else {
    right <- c(left[2:length(left)]-1,n)
  }
  
  BBdata   <- list("data"         = BBobj$data[1:n],
                   "N"            = BBobj$N[1:n],
                   "A"            = BBobj$A[1:n],
                   "left"         = left,
                   "right"        = right,
                   "type"         = type,
                   "params"       = params,
                   "opt"          = BBobj$opt[1:n],
                   "lastchange"   = lastchange+1,
                   "lastparams"   = lastparams,
                   "lasttype"     = lasttype,
                   "pruned"       = BBobj$pruned)
  
  BBobject <- structure(BBdata, class = "BB")
  
  return(BBobject)
}


#Non-homogeneous Poisson Process Generator!!
NHPois <- function(time_length, equation){
  
  expr<-function(x){
    x<-x
    eval(equation)
  }
  current <- 0
  i<-1
  data <- vector()
  maxrate = optim(par=time_length/2,fn=expr,method='L-BFGS-B',lower=0,upper=time_length,control = list(fnscale = -1))$value
  while(current < time_length){
    current<- current+ rexp(n=1,rate=maxrate)
    u <- runif(1)
    
    if (u < expr(current) / ( maxrate)){
      data[i] <- current
      i<-i+1
    }
    
  }
  return(data[c(-length(data))])
}



############################## SCRIPT STARTS HERE ############################################


#create random data
set.seed(1)




equation1 <- expression(50*x+100)
equation2 <- expression(-100*x+250)

equation3 <- expression(25)


x1 <- NHPois(time_length =3,equation = equation1)
x2 <- NHPois(time_length =2,equation = equation2)+x1[length(x1)]
#x3 <- NHPois(time_length =2,equation = equation3)+x2[length(x2)]


x=x1#,x3)#,x4,x5)
#x=x[seq(1,length(x),length.out=length(x)/15)]
#x=x[2:length(x)]
#x=round(x,3)

#x=unique(x)
N <- rep(1,length(x))#let all cells be of size 1

#run model
model <- optinterval(x, N, 5,type=c("power"),pvalue=0.1,verbose=TRUE)

plot(model,show=c("blocks","hist"),bins=100,xlim=c(0,x[length(x)]),ylim=c(0,300),main="Bayesian Blocks")

FREDfit(N,x)


#check summary information and plot
summary(model)


#plot the segmentation in steps
for (i in 1:100){
  png(filename=paste0("~/ne3/",sprintf("%04d",i),".png"),600,500,type='cairo')
  plot(getStep(model,max(which(x<(x[length(x)]*i/100)))),show=c("blocks","points"),binwidth=0.3,ylim=c(0,100),xlim=c(0,x[length(x)]))
  #Sys.sleep(0.1)
  dev.off()
}

dev.off()

save(file="linearmod1.RData",model)
