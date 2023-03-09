# Determine the F critical value

## making genotype x rep a random effect (DFnum = 14; DFden = 14)
# alpha = 0.05
qf(p=.05, df1=14, df2=14, lower.tail=FALSE)
# alpha = 0.001
qf(p=.001, df1=14, df2=14, lower.tail=FALSE)


## using the degrees of freedom from sampling as if they were reps
## Nsample = 5 
N_5<- (5*15*2)-15
N_100<- (100*15*2)-15
# alpha = 0.05
qf(p=.05, df1=14, df2=N_5, lower.tail=FALSE)
qf(p=.05, df1=14, df2=N_100, lower.tail=FALSE)
# alpha = 0.001
qf(p=.001, df1=14, df2=N_5, lower.tail=FALSE)
qf(p=.001, df1=14, df2=N_100, lower.tail=FALSE)



output_all
