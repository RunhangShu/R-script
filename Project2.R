###    Project2  ####
###  Runhang Shu r.shu@ufl.edu ####
###   2019/10/11 ######

# Task 1---------------------------------------
#      Using my name to set myseed and save the dataset  ####
getwd()
setwd("/Users/darson/Desktop/Genomics in R/week8/Week 08")
load(file="dose.data.rdata")
myseed <- sum(grep("[SHU]", LETTERS)) 
set.seed(myseed)
mydata<- dose.data[sample(1:20000,5000,replace = FALSE),]
head(mydata)
save(mydata,file = "Runhang_dose.data.rdata")

###      Load the dataset   ######
load(file = "Runhang_dose.data.rdata")


# Task 2-----------------------------Exploratory analysis-------

#      Conduct PCA analysis by using princomp function 
PCA<-princomp(mydata,cor = TRUE)
head(PCA)
# create my PCA plot function
my_PCA_plot<- function(x, main,col=rep(c(1:20),each=5)) {
    P<-x$loadings
    t.var<-cumsum(x$sdev^2/sum(x$sdev^2))
    var.PC1<-round(t.var[1]*100)
    var.PC2<-round((t.var[2]-t.var[1])*100)
    plot(P[,1],P[,2],
         main = main,
         col=col,
         cex=1, lwd=2,
         xlim=c(min(P[,1])-0.2,max(P[,1])+0.3) ,
         ylim=range(P[,2])*2,
         xlab= paste("PC1:",var.PC1,"% expel.var", sep = ""),
         ylab=paste("PC2:",var.PC2,"% expel.var", sep = "")
        )
    legend(x=0.42,y=0.52,legend=c("Control","100mg","200mg","500mg"),
           col =c(1:4),cex=1,bty="o",pch=1,text.font = 2,text.col = c(1:4),
           pt.cex = 1,pt.lwd = 2,xjust = 0.5,yjust = 0.5)
}
# apply the function I created and input the PCA dataset 
my_PCA_plot(PCA, main="Drugs effect on genes expression")

# positive values on x axis indicate the data is not centered 
mydata.c<-t(apply(mydata, 1, scale, scale=FALSE))
colnames(mydata.c)<-colnames(mydata)
PCA.c<-princomp(mydata.c)
my_PCA_plot(PCA.c, main="Drug doses effect on genes expression of mice")

# find the most responsive genes 
s<-PCA.c$scores
rep.gene<-  s[,1] > 2*PCA.c$sdev[1] & s[,2]< -2*PCA.c$sdev[2]
rep.gene.up<- names(which(rep.gene))


    
    ## Task 3-----------------------------Exploratory analysis-------
# Run t-test between control and 100mg treatment 
c_100<-apply(mydata,1, function(x) {t.test(x[1:5],x[6:10],)$p.value})

# Adjust the p-value by using FDR method 
ad_c_100<-p.adjust(c_100,method = "fdr")
# find the genes that are differentially expressed
DE_100<- names(which(ad_c_100<0.05))
length(DE_100)
## 514 genes were differentially expressed between control-100mg
c_100_higher<-rownames(mydata[apply(mydata[1:5],1,mean)<apply(mydata[6:10],1,mean),])
c_100_UP<- length(which(is.element(c_100_higher,DE_100)))
# 231 genes were upregulated, whilst 283 genes were down-regulated
# Run the same process in control-200mg and control-500mg
c_200<-apply(mydata,1, function(x) {t.test(x[1:5],x[11:15],)$p.value})
ad_c_200<-p.adjust(c_200,method = "fdr")
DE_200<- names(which(ad_c_200<0.05))
length(DE_200)
## 1162 genes were differentially expressed between control-200mg 
c_200_higher<-rownames(mydata[apply(mydata[1:5],1,mean)<apply(mydata[11:15],1,mean),])
c_200_UP<- length(which(is.element(c_200_higher,DE_200)))
# 567 genes were upregualted, whilst 595 genes were downregulated in 200mg treatment 

c_500<-apply(mydata,1, function(x) {t.test(x[1:5],x[16:20],)$p.value})
ad_c_500<-p.adjust(c_500,method = "fdr")
DE_500<- names(which(ad_c_500<0.05))
length(DE_500)
## 1765 genes were differentially expressed between control-200mg 
c_500_higher<-rownames(mydata[apply(mydata[1:5],1,mean)<apply(mydata[16:20],1,mean),])
c_500_UP<- length(which(is.element(c_500_higher,DE_500)))
## 856 genes were upregulated, whilst 909 genes were downregulated in 500mg treatment 


#creat a dataframe of the expression result and draw a barplot 
data.q3<- matrix(c("100mg", "200mg", "500mg",231,567,856,283,595,909),
                 nrow=3,dimnames = list(c(),c("dose","upregulate","downregulate")))

data.q3<-matrix(c(231,283,567,595,856,909),nrow = 2, dimnames=list(c("Upregulate","Downregulate"),c("100mg","200mg","500mg")))
barplot(data.q3,beside = TRUE, horiz = TRUE, main = "Genes Expression Compared with Control",
        xlab = "Number of genes",
        col=c("red","blue"),axis.lty =1,
        legend=rownames(data.q3),
        args.legend = list(x="bottomright"))

### Venn diagram----------------------

# higher in 100mg
high_in_100<-c_100_higher[is.element(c_100_higher,DE_100)]
# higher in 200mg
high_in_200<-c_200_higher[is.element(c_200_higher,DE_200)]
# higher in 500mg
high_in_500<-c_500_higher[is.element(c_500_higher,DE_500)]
# high in both 100mg and 200mg
high_100_200<- high_in_100[is.element(high_in_100,high_in_200)]
# high in both 100mg and 500mg 
hign_100_500<-high_in_100[is.element(high_in_100,high_in_500)]
# high in both 200mg and 500mg 
high_200_500<-high_in_200[is.element(high_in_200,high_in_500)]
# high in all of three treatments 
high_in_all<- high_200_500[is.element(high_200_500,high_in_100)]


#draw a venn plot with upregulated genes in common in treatments 
venn.plot.1 <- draw.triple.venn(
    area1 = 231,
    area2 = 567,
    area3 = 856,
    n12 = 151,
    n23 = 486,
    n13 = 151,
    n123 = 133,
    category = c("100mg", "200mg", "500mg"),
    fill = c("blue", "red", "green"),
    lty = "blank",
    cex = 2,
    cat.cex = 2,
    cat.col = c("blue", "red", "green")
);
grid.draw(venn.plot.1)
grid.newpage()

# calculate the number of downregulated genes in common in treatments 
down_in_100<-DE_100[!(DE_100 %in% high_in_100)]
down_in_200<-DE_200[!(DE_200 %in% high_in_200)]
down_in_500<-DE_500[!(DE_500 %in% high_in_500)]

# down in both 100mg and 200mg
down_100_200<- down_in_100[is.element(down_in_100,down_in_200)] # 164
# down in both 100mg and 500mg 
down_100_500<-down_in_100[is.element(down_in_100,down_in_500)]  # 158
# down in both 200mg and 500mg 
down_200_500<-down_in_200[is.element(down_in_200,down_in_500)]  # 495
# high in all of three treatments 
down_in_all<- down_200_500[is.element(down_200_500,down_in_100)] # 140

#draw a venn plot with downupregulated genes in common in treatments 
venn.plot.2 <- draw.triple.venn(
    area1 = 283,
    area2 = 595,
    area3 = 909,
    n12 = 164,
    n23 = 495,
    n13 = 158,
    n123 = 140,
    category = c("100mg", "200mg", "500mg"),
    fill = c("blue", "red", "green"),
    lty = "blank",
    cex = 2,
    cat.cex = 2,
    cat.col = c("blue", "red", "green"),
    main="ww"
)
grid.draw(venn.plot.2)
grid.newpage()

##### Identify the genes that have a dose-response
as.matrix(mydata)
tapply(q(mydata[1,]), rep(c("control", "100mg", "200mg", "500mg"),each=5), mean)

class(mydata[1,]$Control_2)