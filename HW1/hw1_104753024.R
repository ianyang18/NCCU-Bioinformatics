pam1<-read.table("pam1.txt")

# Amino Acid table's dimension is 20
dim<-dim(pam1)
d<-1:dim

# Normalize each column so that its total sum is equal to 1
for (col in d) {
    pam1[col]<-pam1[col]/sum(pam1[col])
}

# Calculate the PAM25
pam25<-pam1
for (p in 1:25) {
    for (i in d) {
        for (j in d) {
            e<-0
            for (k in d) {
                tmp<-pam25[i,k,]*pam1[k,j,]
                e<-e+tmp
            }
            pam25[i,j]<-e
        }
    }
}

write.table(pam25,"pam25.txt",col.names=TRUE,row.names=TRUE)
