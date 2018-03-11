## Nathan Kraft, UCLA
## nkraft@ucla.edu

## derrived from code by Mark Vellend and Andrew McDonald, from Vellend's Theory of Ecological Communities (Princeton, 2016)
## basic simulation of a two species explicit-metacommunity neutral model with migration but no speciation, adding the potential to introduce fitness differences and frequency dependence that varies among local communties. 


####################################
## main simulation as a function: #
##################################

## first function simulates 2 species in metacommunity, tracking freqency of species 1 over time in output:

simulate_metacommunity=function(years=50, patches=10, J=100, m=0, fitness_ratio_ave=rep(1,patches), frequency_dependence=rep(0,patches)){
	
	###############################
	## Description of Arguments: #
	#############################
	
	## years- here definied arbitrarily as (J * patches) deaths = 1 year
	## this is helpful if you want to vary J and compare results on same timescale
	
	## patches- number of local communities to track
	
	## J- number of individuals per patch. 
	## Here the metacommunity is modeled explicitly, so Jm (metacommunity size) is:
	Jm <- J * patches
	
	## m- migration rate, probability that a local birth
	## comes from the metacommunity instead of a local individual
	## needs to be between 0 and 1. 
	
	
	## fitness_ratio_ave specifies the fitness difference between
	## species 1 and 2 in each patch.
	## Needs to be vector of length 'patches'
	## GUIDANCE ON VALUES:
	## 1= no fitness difference between species 1 and 2 (neutral)
	## if you want symmetric fitness differences in different patches
	## your values need to be INVERSES of eachother, such as c(1.2, 1/1.2), NOT c(1.2, 0.8)
	
	
	## frequency_dependence specifies the slope of frequency dependence
	## in each patch. Needs to be vector of length 'patches"
	## 0= no frequency dependence. 
	## Note that fitness_ratio_average above determines the intercept of the 
	## relationship for frequency dependence 
	
	### a few failsale tests to stop the simulation if key arguments are outside of spec:
	if(m<0 | m>1){return(print("simulation canceled: m needs to be between 0 and 1"))}
	if(length(fitness_ratio_ave)!=patches){return(print("simulation canceled:fitness_ratio_ave needs to be length 'patches' "))}
	if(length(frequency_dependence)!=patches){return(print("simulation canceled:fitness_ratio_ave needs to be length 'patches' "))}
	
	########################
	### start simulation ##
	######################
	
	## initialize output matrix, where rows are years and columns are communities (patches). 
	##We'll be working with two species in a zero sum model, so we only need to track one species
	frequency_sp1 <- matrix(nrow = years, ncol = patches) 
	
	## initial population size of species 1- set to be equally abundand with species 2 to start
	initial_abundance_sp1 <- 0.5*J 
	
	## our live, working snapshot of the metacommunity, where rows are individuals and columns are local communities. 
	metacommunity_state <- matrix(nrow=J, ncol=patches)
	
	## set half the individuals in all patches to species 1:
	metacommunity_state[1:initial_abundance_sp1,] <- 1
	## and the other half to species 2:
	metacommunity_state[(initial_abundance_sp1+1):J,] <- 2 
	
	
	## record frequecy of species 1 for year 1- same for all patches initially.
	frequency_sp1[1,]<-initial_abundance_sp1/J
	
	## for clarity here I will nest 'for' loops here for years and deaths, though this makes it slower to run:
	
	for(i in 2:years){
		## assuming our intial data represents year 1, we'll start at year 2:
		
		for(j in 1:Jm){
			## each year has J * patches (Jm) deaths in it- we'll loop over that here:
			
			## for each death, choose one patch as location- all have same chance:
			patch<-sample(1:patches, 1)
			
			##calculate probabilty of replacement by species 1:
			
			if(runif(1)<m){
				## randomly choose a number between 0 and 1,
				## if it is less than m, death is replaced by
				## dispersal from metacommunity.
				## chance of replacement being species 1 is just 
				## the relative abundance of species 1 in species pool
				Prob_repro_Sp1<-sum(metacommunity_state==1)/Jm
			} else {
				## if not, replacement is from local community.
				## fitness differences and frequency dependence, if present
				## play out in the local replacement process here
				
				# calculate frequencies of species 1 and 2 in selected patch:
				local_freq_1 <- sum(metacommunity_state[,patch]==1)/J
				local_freq_2 <- 1- local_freq_1
				
				## caculates the fitness ratio between sp 1 and 2 in the randomly 
				## selected patch based on the initial vector of frequency dependences
				## and the fitness ratio averages in each patch
				
				fitness_ratio <- exp(frequency_dependence[patch]*(local_freq_1-0.5)+log(fitness_ratio_ave[patch]))
				
				#  use this ratio and oberved frequencies of species to calculate chance that species 1 reproduces:
				Prob_repro_Sp1 <- fitness_ratio*local_freq_1 / ( fitness_ratio*local_freq_1 + local_freq_2)
				
			}
			
			## now select a random individual inside the selected patch to be killed
			## replacement is by sp 1 at probablility calculated above via sample()
			metacommunity_state[ ceiling(J*runif(1)) , patch] <-sample(c(1,2),1, prob=c(Prob_repro_Sp1, 1- Prob_repro_Sp1))
			
			
		}
		
		## at the end of each "year" (= Jm deaths), record state of metacommunity:
		frequency_sp1[i,] <-colSums(metacommunity_state==1)/J ## counts frequency of sp 1 in each patch
		
		
	}
	
	## after the years are done, return the final dataset on frequency of sp 1 over time in each patch:
	return(frequency_sp1)

}



#########################################################################
## function for simple line plots of output from above simulation: ###
####################################################################


trace_metacommunity=function(df, title="your title here"){
	## expects output from simulate_metacommunity()
	## you can add any title you like to track parameters, etc
	plot(1:nrow(df), df[,1], type='l', xlab="Years", ylab="Freq. of Species 1", ylim=c(0,1), main=title)
	for(i in 2:(ncol(df))){
		lines(1:nrow(df), df[,i], type='l')
	}
}


###################################
## simple examples ####
####################

## run functions above first to save into memory.

#pdf("Sims_basic processes.pdf", width=8, height=9)
#layout(matrix(1:6, nrow=3, byrow=T))

### local processes: drift with no migration:
#simulate_metacommunity(years=50, patches=10, J=100, m=0, fitness_ratio_ave=rep(1,10), frequency_dependence=rep(0,10))->output
#trace_metacommunity(output, "neutral, J=100, m=0")

## effect of J on drift, still no migration:
#simulate_metacommunity(years=50, patches=10, J=1000, m=0, fitness_ratio_ave=rep(1,10), frequency_dependence=rep(0,10))->output
#trace_metacommunity(output, "neutral, J=1000, m=0")

## constant selection- e.g. habitat filtering:

## one patch type, no migration:
#simulate_metacommunity(years=50, patches=10, J=100, m=0, fitness_ratio_ave=rep(1.2,10), frequency_dependence=rep(0,10))->output
#trace_metacommunity(output, "constant selection, fitRatio = 1.2, J=100, m=0")

## two symmetric patch types, no migration:
#simulate_metacommunity(years=50, patches=10, J=100, m=0, fitness_ratio_ave=c(rep(1.2,5), rep(1/1.2, 5)), frequency_dependence=rep(0,10))->output
#trace_metacommunity(output, "constant selection, 2 habitats, fitRatios = 1.2 & 1/1.2")

## weak negative frequency dependence with a fitness difference, no migration:

#simulate_metacommunity(years=50, patches=10, J=100, m=0, fitness_ratio_ave=rep(1.1,10), frequency_dependence=rep(-.3,10))->output
#trace_metacommunity(output, "freq dependence -.3, fitRatio = 1.1, J=100, m=0")

## strong negative frequency dependence with a fitness difference, no migration:

#simulate_metacommunity(years=50, patches=10, J=100, m=0, fitness_ratio_ave=rep(1.1,10), frequency_dependence=rep(-.8,10))->output
#trace_metacommunity(output, "freq dependence -.8, fitRatio = 1.1, J=100, m=0")

#dev.off()

