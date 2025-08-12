library(tibble)
# Differential proportion analysis of cell populations between conditions
source("D:/program/DPA/diffprop_functions.R")
normal=subset(S1, subset = group =="dbm")
dkd=subset(S1, subset = group =="dbdb")
a=as.data.frame(table(normal@active.ident))
b=as.data.frame(table(dkd@active.ident))
colnames(a)=c("cellstate","CTL")
colnames(b)=c("cellstate","DKD")
data=merge(a,b,by="cellstate")
data=column_to_rownames(data, var = "cellstate") %>% t()
res.table = c()

## Go through a series of error probabilities
for (err_prob in c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05)) {
  tip.exp <- generateNull(data, n=100000, p=err_prob);
  res.1 = two.class.test(data, tip.exp, cond.control="CTL", cond.treatment="DKD",to.plot=F)
  res.table = rbind(res.table, res.1)}
rownames(res.table) = as.character(c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05))
write.csv(res.table, file="D:/DPA.csv")
