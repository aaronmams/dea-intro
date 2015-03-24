#====================================================================================================
#Solve the output oriented DEA problem


#The function dea.out will use Rglpk to solve the output oriented DEA program.  The function 
# returns a technical efficiency measure between 0 and 1 reflecting the level of technical
# efficiency for each DMU. A score of 1 indicates technical efficiency.

# Function takes the following inputs:
#  data 			--> A data frame, should be organized such that each DMU is a row
#					and columns include inputs and outputs  	
#	output.names  	--> A character vector including the column names of data indicating
#					the output columns
#	input.names 	--> A character vector includinig the column names of data indicating 
#					input columns
#	rts			--> A string declaring returns to scale.  rts = 'crs' solves the constant
#					returns to scale problem.  rts = 'vrs' solves the variable 
#					returns to scale problem.

# Function returns the following outputs:
#
#	te 			--> Matrix.  1st column is the technical efficiency score. 2nd column is the 'status'.  From
#					the Rglpk documentation: as long as the control parameter is set to the default
#					status returns "0" if an optimal solution is found and non-zero integer otherwise.

library(Rglpk)

dea.out <- function(data,output.names, input.names, rts, slack){
  ndmu <- nrow(data)
  obj <- c(1,rep(0,ndmu))
  
  s <- length(output.names)
  m <- length(input.names)
  z <- matrix(0,nr=ndmu,nc=ndmu+s+m+3)
  
#======================================================================================================
#======================================================================================================
#======================================================================================================
#======================================================================================================

  #set up the output matrix for the function...te will store the solution to the LP for each DMU
  # and an integer reflecting the status of the optimization
  te <- matrix(0,nr=ndmu,nc=ndmu+2)
  
  for(irow in 1:nrow(data)){
    #set up y matrix
    y.dmu <- as.matrix(data[irow,which(colnames(data)%in% output.names)])
    yrest <- data[1:ndmu,which(colnames(data)%in% output.names)]
    y <- cbind(t(-y.dmu),t(yrest))
    
    #set up x matrix
    x.dmu <- as.matrix(data[irow,which(colnames(data) %in% input.names)])
    xrest <- data[1:ndmu, which(colnames(data) %in% input.names)]
    x <- cbind(matrix(0,nr=length(input.names),nc=1),t(xrest))
    
    #set up nonzero constraints on the lambdas
    lambda.mat <- cbind(matrix(0,nr=ndmu,nc=1),diag(ndmu))
    
    #Max or Min
    max = TRUE
    
    #Set up direction of inequality & RHS.
    #variable returns to scale means the sum of lambdas = 1
    #constant returns to scale just imposes the non-zero constraint on each of the lambdas
    if(rts == "vrs"){
      #combine them into the coefficient matrix
      A <- rbind(y,x,lambda.mat,cbind(0,matrix(1,nr=1,nc=ndmu)))
      rhs <- c(rep(0,ncol(y.dmu)),t(x.dmu),rep(0,ndmu),1)
      dir <- c(rep('>=',nrow(y)),rep('<=',nrow(x)),rep('>=',nrow(lambda.mat)),"==")
      result <- Rglpk_solve_LP(obj=obj,mat=A,dir=dir,rhs=rhs,types=c(rep("C",ndmu+1)),max=max)
      
    }
    if(rts == "crs"){
      #combine them into the coefficient matrix
      A <- rbind(y,x,lambda.mat)
      rhs <- c(rep(0,ncol(y.dmu)),t(x.dmu),rep(0,ndmu))
      dir <- c(rep('>=',nrow(y)),rep('<=',nrow(x)),rep('>=',nrow(lambda.mat)))
      result <- Rglpk_solve_LP(obj=obj,mat=A,dir=dir,rhs=rhs,types=c(rep("C",ndmu+1)),max=max)
      
    }
    
    te[irow,] <- c(result$solution, result$status)
    
  }
  te <- as.data.frame(te)
  names(te) <- c("theta",paste("lambda",1:ndmu,sep=""),"optim.code")

#======================================================================================================
#======================================================================================================
#======================================================================================================
#======================================================================================================
# if SLACK = T then solve the second stage problem

if(slack==T){
  
  #-----------------------------------------------------------------------------------------------------------------------
  #2nd Stage: First project the observed output bundles to the frontier:
  y.new <- data.frame(data[,which(colnames(data) %in% output.names)] * te[,which(names(te)=="theta")])
  
  #Next maximize the sum of slacks given the efficient output bundle:
  obj <- c(rep(0,ndmu),rep(1,s+m))
  
  #attach slack parameters to output rows
  A <- cbind(t(data[,which(colnames(data) %in% output.names)]), -diag(s), matrix(0,nr=length(output.names),nc=m))
  
  #attach slack parmeters to input rows and rbind with output rows
  A <- rbind(A, cbind(t(data[,which(colnames(data) %in% input.names)]),matrix(0,nr=length(input.names),nc=s),diag(m)))
  
  for(idmu in 1:ndmu){
    y.now <- y.new[idmu,]
    x.now <- data[idmu,which(colnames(data) %in% input.names)]
    
    if(rts=='crs'){
      AA <- rbind(A, cbind(diag(ndmu),matrix(0,nr=ndmu,nc=s+m)))
      #establish direction of constraint
      dir <- c(rep('==',s+m),rep('>=',ndmu))
      rhs <- c(unlist(c(y.now,x.now)),rep(0,ndmu))
      
    }
    
    if(rts=='vrs'){
      #establish equality constraint for lambdas and lower bounds for lambdas and slacks
      AA <- rbind(A, cbind(matrix(1,nr=1,nc=ndmu),matrix(0,nr=1,nc=s+m)))
      #establish direction of constraint
      dir <- c(rep('==',s+m),'==')
      rhs <- c(unlist(c(y.now,x.now)),1)
      
    }
    
    
    z.tmp <- Rglpk_solve_LP(obj=obj,mat=AA,dir=dir,rhs=rhs,types=c(rep('C',ndmu+s+m)),max=T)
    z[idmu,1] <- theta[idmu,1]
    z[idmu,2:ncol(z)] <- c(z.tmp$solution,z.tmp$optimum, z.tmp$status)
  }
  #---------------------------------------------------------------------------------------------------------------
  #Save output and format a little:
  
  z <- as.data.frame(z)
  names(z) <- c("theta",paste("lambda",1:ndmu),rep("Sy",s),rep("Sx",m),"optimum","optim.code")
  return(z)
}else{
  return(te)
}


#======================================================================================================
#======================================================================================================
#======================================================================================================
#======================================================================================================

}
#=====================================================================================================


