# Grace Chung
# Tau-Leaping SCIR model 

# Setup Times 
tMax = 1460               # Total time to simulate
tau = 0.1                 # Time step to jump by
t = seq(0,tMax,tau)       # Array of times we will simulate

# Parameters 
N = 40000;

beta = 0.00085;
contact = 10.86;
gammac = 0.00658;
delta = 0.000870;
mu = 0.00000000315;
gamma = 0.152;
recovery = 0.00137;

# NO VACCINATION
# Initial Conditions
C0 = 1;
I0 = 1;
S0 = N-I0-C0;
R0 = 0;

# Array to store model simulation results
x = matrix(0,length(t),4); 
x[1,] = c(S0, C0, I0, R0); 

incidence = matrix(0,length(t),1)

y = matrix(0,1000,1)

# We will suppose the number of events within each timestep tau is Poisson-distributed with
# mean given by the total number of events expected per timestep from the above equations. 
for (j in 1:1000){
  for (i in 1:(length(t)-1)){
    NewCarrier = rpois(1, tau*(beta*contact/N)*(x[i,1])*(0.5*x[i,2]+x[i,3]));
    CarriersRecovered = rpois(1, tau*gammac*x[i,2]);
    CarriersProgress = rpois(1, tau*delta*x[i,2]);
    InfectedDie = rpois(1, tau*mu*x[i,3]);
    CurrRecov = rpois(1, tau*gamma*x[i,3]);
    Recov = rpois(1, tau*recovery*x[i,4]);
    
    if(NewCarrier > x[i,1]){  #if number of new carriers > number currently susceptible
      NewCarrier = x[i,1];    #then only the max possible could have switched before the variable hits zero
    }
    
    if((CarriersProgress + CarriersRecovered) > x[i,2]){  
      CarriersProgress = (delta/(delta+gammac))*x[i,2];    
      CarriersRecovered = (gammac/(delta+gammac))*x[i,2];
    }
    
    if((InfectedDie + CurrRecov) > x[i,3]){  
      InfectedDie = (mu/(mu+gamma))*x[i,3];    
      CurrRecov = (gamma/(mu+gamma))*x[i,3];
    }
    
    if(Recov > x[i,4]){  #if number of new carriers > number currently susceptible
      Recov = x[i,4];    #then only the max possible could have switched before the variable hits zero
    }
   
    x[i+1,] = x[i,] + c(-NewCarrier+CarriersRecovered+Recov, NewCarrier-CarriersProgress-CarriersRecovered, CarriersProgress-InfectedDie-CurrRecov, CurrRecov-Recov);
    
    incidence[i,] = c(CarriersProgress);
  } 
  
  y[j,] = c(sum(incidence[1:14600,1]))
}  
  
sum(y[,1])
sum(y[,1])/160000000

write.csv(y, "no vaccination.csv")


# VACCINATION
# Initial Conditions
p = 0.90
q = 0.73

C0 = 1;
I0 = 1;
S0 = (N-C0-I0)*(1-(p*q));
R0 = 0;

# Array to store model simulation results
x = matrix(0,length(t),4); 
x[1,] = c(S0, C0, I0, R0); 

incidence = matrix(0,length(t),1)

y = matrix(0,1000,1)

# We will suppose the number of events within each timestep tau is Poisson-distributed with
# mean given by the total number of events expected per timestep from the above equations. 
for (j in 1:1000){
  for (i in 1:(length(t)-1)){
    NewCarrier = rpois(1, tau*(beta*contact/N)*(x[i,1])*(0.5*x[i,2]+x[i,3]));
    CarriersRecovered = rpois(1, tau*gammac*x[i,2]);
    CarriersProgress = rpois(1, tau*delta*x[i,2]);
    InfectedDie = rpois(1, tau*mu*x[i,3]);
    CurrRecov = rpois(1, tau*gamma*x[i,3]);
   
    if(NewCarrier > x[i,1]){  #if number of new carriers > number currently susceptible
      NewCarrier = x[i,1];    #then only the max possible could have switched before the variable hits zero
    }
    
    if((CarriersProgress + CarriersRecovered) > x[i,2]){  
      CarriersProgress = (delta/(delta+gammac))*x[i,2];    
      CarriersRecovered = (gammac/(delta+gammac))*x[i,2];
    }
    
    if((InfectedDie + CurrRecov) > x[i,3]){  
      InfectedDie = (mu/(mu+gamma))*x[i,3];    
      CurrRecov = (gamma/(mu+gamma))*x[i,3];
    }
    
    x[i+1,] = x[i,] + c(-NewCarrier+CarriersRecovered, NewCarrier-CarriersProgress-CarriersRecovered, CarriersProgress-InfectedDie-CurrRecov, CurrRecov);
    
    incidence[i,] = c(CarriersProgress);
  } 
  
  y[j,] = c(sum(incidence[1:14600,1]))
}    
  
sum(y[,1])
sum(y[,1])/160000000

write.csv(y, "90% vaccination 73% vaccine efficacy.csv")