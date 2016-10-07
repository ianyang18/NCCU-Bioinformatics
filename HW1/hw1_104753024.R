##########################################
# Goal: Calculating the PAM250 matrix.
# Input data: pam1.txt
# Output file: pam250.txt
##########################################

# library "gdata" is suitable for data alignment, adjusting the data formation
if(!require("gdata")) {
    install.packages("gdata",repos="http://cran.rstudio.com/")
    library("gdata")
} else {
    library("gdata")
}

# Read the data from "pam1.txt"
pam1<-read.table("pam1.txt")

# Amino Acid table's dimension is 20 stored in dimension
dimension<-dim(pam1)
d<-1:dimension

# Normalize each column so that its total sum is equal to 1
for (col in d) {
    pam1[col]<-pam1[col]/sum(pam1[col])
}

# Calculate the PAM250 (Matrix Multiplication)
# Variable:
#   pam250 is a data.frame preserving the final result
#   p is the variable counting the multiplication times.
#   i, j is the row, column index of the PAM matrix.
#   e is the element storing the middle result of the inner product.
pam250<-pam1
for (p in 1:250) {
    for (i in d) {
        for (j in d) {
            e<-0
            for (k in d) {
                tmp<-pam250[i,k,]*pam1[k,j,]
                e<-e+tmp
            }
            pam250[i,j]<-e
        }
    }
}

# Round the number to decimal places 2
pam250r<-round(pam250, digits = 2)
pam250s<-pam250r*100

# Write the result into the file "pam250.txt"
write.table(pam250,"pam250.txt",col.names=TRUE,row.names=TRUE)
write.table(pam250r,"pam250_round.txt",col.names=TRUE,row.names=TRUE)
write.fwf(pam250s,file="pam250_scale100.txt",sep="\t",quote=F,rownames=T,rowCol=" ",justify="right")
