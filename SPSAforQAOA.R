library("QuantumOps")

#Create QAOA circuit (represented by a matrix) from edges and parameters
QAOA <- function(nQubits,edges,p,gammas,betas){
    #Initial Hadamard gates
    m <- many(gate=H(),n=nQubits)
    #For each stage
    for(P in 1:p){
        #Phase separation for each edge
        for(e in 1:dim(edges)[1]){
            m <- controlled(gate=X(),n=nQubits,cQubits=edges[e,1],tQubit=edges[e,2]) %*% m   #CNOT from control to target
            m <- single(gate=Rz(gammas[P]),n=nQubits,t=edges[e,2]) %*% m                     #Rz on target
            m <- controlled(gate=X(),n=nQubits,cQubits=edges[e,1],tQubit=edges[e,2]) %*% m   #CNOT from control to target
        }
        #Rx on all qubits
        m <- many(gate=Rx(2*betas[P]),n=nQubits) %*% m
    }
    m #Return the matrix
}

#classically check the score (cut) of a bitstring
computeScore <- function(bitstring,edges){
    score <- 0
    for(e in 1:dim(edges)[1]){  #For each edge  (R indices start at 1)
        if(bitstring[edges[e,1]+1] != bitstring[edges[e,2]+1]) #If bits on either end of the edge are different
            score <- score + 1                             #Score increases by 1
    }
    score #return score
}

#Find the expectation value of QAOA for a given set of parameters
ExpectationValue <- function(n,edges,p,gammas,betas,nSamples){
    #Generate the QAOA quantum circuit 
    circuit <- QAOA(nQubits=n,edges=edges,p=p,gammas=gammas,betas=betas)
    #Initialze a qubit register in the |00...0> state
    qr <- intket(x=0,n=n)
    #Apply the circuit (represented by a matrix) to the quantum register (represented by a vector)
    qr <- circuit %*% qr
    #Simulate the measurements of the quantum state (each measurement would correspond to individual run in experiment)
    measurements <- sample(x=0:(2^n-1),size=nSamples,prob=probs(qr),replace=TRUE) #Sample from 0 to n-1 with probabilities from the quantum state
    
    #Test the score of the measurements
    scoreSum <- 0     #Add up all the scores from all measurments
    for(i in 1:nSamples)  #convert integer result to binary bitstring before parsing
        scoreSum <- scoreSum + computeScore(bitstring=convert_dec2bin(measurements[i],len=n),edges=edges)
    expectationValue <- scoreSum/nSamples #Average them
    
    return(expectationValue)
}

#Run Experiment
SPSAforQAOA <- function(n,edges,p,nIterations,nSamples,a_start,c_start,decay=2){
    #Make learning rate and perturbation sequences
    a <- a_start / (1:nIterations)^decay
    c <- c_start / (1:nIterations)^decay
    c[c<0.01] <- 0.01 #put a lower bound on c (don't make perturbations too small)
    #Generate initial parameters at random
    gammas <- runif(min=-0.05,max=0.05,n=p)
    betas  <- runif(min=-0.05,max=0.05,n=p)
    
    #Run iterations of SPSA
    for(i in 1:nIterations){
        #Generate perturbation vectors (bournoulli random variables of magnitude c[i])
        dGammas <- sample(x=c(-1,1),size=p,replace=TRUE)*c[i]
        dBetas  <- sample(x=c(-1,1),size=p,replace=TRUE)*c[i]
        
        #Compute the expectation values of the + and - configurations
        Fplus  <- ExpectationValue(n=n,edges=edges,p=p,gammas=gammas+dGammas,betas=betas+dBetas,nSamples=nSamples)
        Fminus <- ExpectationValue(n=n,edges=edges,p=p,gammas=gammas-dGammas,betas=betas-dBetas,nSamples=nSamples)
        
        #Compute the gradients of the parameters based on expectation values
        gGammas <- (Fplus - Fminus) / (2*dGammas)
        gBetas  <- (Fplus - Fminus) / (2*dBetas)
        
        #Update the parameters based on gradients
        gammas <- gammas + a[i] * gGammas
        betas  <- betas  + a[i] * gBetas
        
        #Report the progress
        print(paste('Iteration:',i,'Exp(+):',sprintf('%.4f',Fplus),'Exp(-):',sprintf('%.4f',Fminus)))
    }
}

#The provided example problem
example <- function(){
    #4 qubits
    n <- 4
    #Edges of the max-cut problem in matrix from (rows are edges, column 1 is source vertex, column 2 is destination vertex)
    edges <- matrix(c( 0,1, 1,2 , 2,3 , 0,3),byrow=TRUE,ncol=2)
    #p=2 is sufficient for this problem
    p <- 2
    #A sufficient number of optimization iterations to solve problem
    nIterations <- 100
    #Typically need quite a few samples (measurements of quantum circuit) per iteration to 
    nSamples <- 10000
    #Heuristically chosen a and c
    a_start <- 0.25
    c_start <- 0.25
    decay <- 0.5
    SPSAforQAOA(n=n,edges=edges,p=p,nIteration=nIterations,nSamples=nSamples,a_start=a_start,c_start=c_start,decay=decay)
}