# retry the fourier expansion for the encode process was wrong!

rm(list = ls())
#library(BiocInstaller)
#biocLite('FRGEpistasis')

#library(FRGEpistasis)

library(fda)
library(MASS)
library(stats)
library(factoextra)


#=================function define by ourself======================================
#new defined position
fourierExpansion <-function(gene_idx,geno,gene_list,snp_map,rng){
  #missing(geno)
  snpName<-as.vector( snp_map[,2] )
  snp.position<-as.vector( snp_map[,4] )
  snp.chrom=snp_map[,1]
  
  #gene.chrom<-as.vector( gene_list$Chromosome )
  gene.symbol<-as.vector( gene_list$Gene_Symbol)
  
  #gene.start<-subset(gene_list,Gene_Symbol==gene_idx)$Start
  #gene.end<-subset(gene_list,Gene_Symbol==gene_idx)$End
  #gene.end<-as.vector(gene_list$End
  geno<-as.data.frame(geno)
  
  #sub.idx<-as.character(subset(newposi,Genenames==gene_idx)$SNPnames)
  #glist<-gene_list$Gene_Symbol
  sub.idx<-glist%in%gene_idx
  sub.idx<-names(geno)[sub.idx]
  #sub.idx<-as.numeric(sub.idx)
  #sub.idx<-snp.chrom == gene_idx
  #sub.idx<-sub.idx & (snp.position>(gene.start -rng) )
  #sub.idx<-sub.idx & (snp.position<(gene.end+rng) )
  #scope<-(snp.position>(gene.start -rng) )&&(snp.position<(gene.end+rng) )
  #if(sum(sub.idx)>0){
  #x<-as.matrix(geno[,sub.idx])
  #  print('get the right locis')
  #}else{sub.idx<-snp.chrom == gene_idx
  #}
  x<-as.matrix(geno[,sub.idx])
  #pos<-as.character(subset(newposi,Genenames==gene_idx)$Position)
  #posivec<-as.vector(newposi$Position)
  pos<-posivec[sub.idx]
  pos<-as.numeric(pos)
  x<-as.data.frame(x)
  for(i in 1:dim(x)[2]){
    x[,i]<-as.numeric(as.character(x[,i]))
  }
  x<-as.matrix(x)
  
  nsnps <-dim(x)[2]
  
  gil=list();
  # if the number of snp in the gene is over 3, expansion occurs
  if(nsnps>1)
  {
    if ( is.null(pos) )
    {
      pos <-(0:( nsnps-1) )/(nsnps-1)
    }else {
      idx<-order(pos)
      print(pos)
      print(idx)
      print(colnames(x))
      #if(length(order)!=1){
      x<-x[,idx]
      
      pos<-pos[idx]
      pos<-(pos-pos[1])/(pos[nsnps]-pos[1])
    }
    
    # use PCA to determine the number of basis
    eigenval<-prcomp(x)$sd^2
    sum_eigen=sum(eigenval)
    tmp=0
    n_of_basis=0
    for(i in 1:length(eigenval))
    {
      tmp=eigenval[i]+tmp
      n_of_basis=i;
      if(tmp>=0.8*sum_eigen) break
    }
    n_of_basis=floor(n_of_basis/2)*2+1 
    #make n_of_basis_A the odd number
    # end of setting the basis number
    
    frange <-c(pos[1], pos[length(pos)])
    
    fbasis<-create.fourier.basis(frange,nbasis=n_of_basis)
    phi=eval.basis(pos,fbasis);
    gil$zeta <-t(ginv(t(phi)%*%phi)%*%t(phi)%*%t(x))
  }else gil$zeta  <-x
  
  return(gil$zeta)
  
}

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#########==================below is the demostration of FRGEpistasis

##read the list of genotype files
dir_path = '/Users/ajmd/code/R/newFourier/data'
geno_info <- read.table('/Users/ajmd/code/R/newFourier/data/simGeno-chr2.raw',header=T)

##read the gene annotation file
#gLst<-read.csv(system.file("extdata", "gene.list.csv", package="FRGEpistasis"))
gLst = read.csv(paste(dir_path, 'gene.list.csv', sep = '/'))

##define the extension scope of gene region
#rng=0
#fdr=0.05

## output data structure
#out_epi <- data.frame( )
##log transformation
#phenoInfo [,2]=log(phenoInfo [,2])
##rank transformation
#out_epi = fRGEpistasis(work_dir,phenoInfo,geno_files,mapFiles,gLst,fdr,rng)
#c=0.5
#phenoInfo[,2]=rankTransPheno(phenoInfo[,2],c)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



#import new data
# modify the path to source your data
library(snpStats)
#frpdir<-"C:/Users/jiaxin/Documents/R/win-library/3.3/GeneGeneInteR/extdata/"
#ped<-"C:/Users/jiaxin/Documents/R/win-library/3.3/GeneGeneInteR/extdata/inter.ped"
ped = paste(dir_path,'inter.ped',sep = '/')
snpinfo<-read.table(ped,header = F)
#info <-"C:/Users/jiaxin/Documents/R/win-library/3.3/GeneGeneInteR/extdata/inter.info"
info = paste(dir_path, 'inter.info',sep = '/')
newped<-snpStats::read.pedfile(file=ped, snps = info)

names(snpinfo)[1:6]<-names(geno_info)[1:6]


#posi <-"C:/Users/jiaxin/Documents/R/win-library/3.3/GeneGeneInteR/extdata/inter.txt"
posi = paste(dir_path, 'inter.txt',sep = '/')

newposi<-read.table(posi,header = T)
newmap<-newposi
newmap[,1]<-newposi[,3]
newmap[,3]<-rep(0,dim(newmap)[1])
names(newmap)<-c('chr','snp_id','distant','posi')
#so the compute map file we can get from below
newmap<-newmap[order(newmap$chr),]
newmap$snp_id<-as.character(newmap$snp_id)
for(i in 1:dim(newmap)[1]){
  newmap[i,2]<-gsub('rs','r',newmap[i,2])
}
#write.csv(newmap,'c:/Users/jiaxin/Desktop/??????ϵ????/sorteddata/newmap.csv',row.names = F)


#we want get the snp info
{
  #casesnp<-read.csv('c:/Users/jiaxin/Desktop/??????ϵ????/newfrge/case.csv')
  casesnp = read.csv(paste(dir_path, 'case.csv', sep = '/'))
  casesnp$X.3<-rep(1,dim(casesnp)[1])
  #ctrlsnp<-read.csv('c:/Users/jiaxin/Desktop/??????ϵ????/newfrge/control.csv')
  ctrlsnp = read.csv(paste(dir_path, 'control.csv', sep = '/'))
  ctrlsnp$X.3<-rep(0,dim(ctrlsnp)[1])
  
  names(ctrlsnp)[6:143]<-gsub('X','rs',names(ctrlsnp)[6:143])
  names(casesnp)[6:143]<-gsub('r','rs',names(casesnp)[6:143])
  snpmatrix<-rbind(casesnp,ctrlsnp)
  newcol<-as.data.frame(matrix(rep(0,dim(snpmatrix)[1]),ncol = 1))
  rownames(newcol)<-rownames(snpmatrix)
  #names(newcol)<-names(geno_info)[5]
  fisnp<-cbind(snpmatrix[,1:4],newcol,snpmatrix[,-c(1:4)])
  names(fisnp)[1:6]<-names(geno_info)[1:6]
  
  #recode the snp with minor allele frequence
  #load the prepared function lociencode
  #source('c:/Users/jiaxin/Desktop/??????ϵ????/newcasuality/recode_minor_Freq.R')
  source(paste(dir_path, 'recode_minor_Freq.R', sep = '/'))
  
  
  snpm<-lociencode(fisnp[,-c(1:6)])
  snpm<-cbind(fisnp[,1:6],snpm)
  }
snpm$IID<-snpm$FID
#write.csv(snpm,'c:/Users/jiaxin/Desktop/??????ϵ????/newfrge/snpm.csv',row.names = F)



#also phenoinfo is needed
newpheno<-snpm[,c(2,6)]

#then we need the annotation file gLst
genelist<-as.data.frame(table(newposi$Genenames))
#genebed<-read.csv('c:/Users/jiaxin/Desktop/??????ϵ????/newfrge/mapfile.csv')
genebed = read.csv(paste(dir_path, 'mapfile.csv',sep = '/'))

newglst<-cbind(genelist,genebed)[,-2]
names(newglst)<-names(gLst)

#newglst<-newglst[order(newglst$Chromosome),]
#save 
#write.csv(newglst,'c:/Users/jiaxin/Desktop/??????ϵ????/sorteddata/newglst.csv',row.names = F)


##define the extension scope of gene region
rng=0
fdr=0.05
## output data structure
out_epi <- data.frame( )


# 
##log transformation
#newpheno[,2]=log(newpheno[,2])
##rank transformation
#out_epi = fRGEpistasis(work_dir,phenoInfo,geno_files,mapFiles,gLst,fdr,rng)
#c=0.5
#phenoInfo[,2]=rankTransPheno(phenoInfo[,2],c)
#library(FRGEpistasis)
#genefre<-fRGEpistasis(wDir=frpdir,phenoInfo=newpheno,gnoFiles=snpm,
                     # mapFiles=newmap,gLst=newglst,fdr,rng)



#the fourier expansion performed below
#gene_idx<-gIdx_i;geno<-gnoChr;gene_list<-gListChr;snp_map<-mapChr#rng
#newpheno<-read.csv('c:/Users/jiaxin/Desktop/??????ϵ????/newfrge/')

#load the gene and position information
#gene_info<-read.csv('c:/Users/jiaxin/Desktop/??????ϵ????/sorteddata/gene_info.csv')
gene_info = read.csv(paste(dir_path, 'gene_info.csv',sep = '/'))

#formercase<-read.csv('c:/Users/jiaxin/Desktop/??????ϵ????/sorteddata/finalsnpm.csv')
formercase = read.csv(paste(dir_path, 'finalsnpm.csv', sep = '/'))

#formercase<-cbind(pheno$PHENOTYPE,formercase);names(formercase)[1]<-'pheno'
missing<-names(formercase)%in%gene_info$SNPnames
miss_snp<-names(formercase)[!missing]
totalsnp<-c(as.vector(gene_info$SNPnames),miss_snp)
#totalsnp<-c(as.vector(gene_info$SNPnames),miss_snp)
searchedname<-'MAOA MAOA SLC1A2 NR3C1 CYP2D6 LOC100288866 GRIN2A AKT3 BDNF FAM MAPK8 GRM7 GAD1 GRIN2B ABCB6 GRIK4 ABCB1 NTRK2 ABCB1 NULL NULL MC2R ABCB1 CRHBP CYP2C9'
searchedname<-strsplit(searchedname,split = ' ')
totalgene<-c(as.character(gene_info$Genenames),searchedname[[1]])
#gene_info$SNPnames<-totalsnp
rawvec<-totalgene
names(rawvec)<-totalsnp

newgeno<-snpm[,-c(1:6)]
glist<-as.character(rawvec[names(newgeno)])
#newposi<-subset(newposi,SNPnames%in%names(newgeno))

#define a posi vecector instead use the newposi
posivec<-as.vector(newposi$Position)
names(posivec)<-as.vector(newposi$SNPnames)



geno_expn=matrix(0,dim(snpm)[1],1)
geno_expn=geno_expn[,-1]
#newglst<-read.csv('c:/Users/jiaxin/Desktop/??????ϵ????/sorteddata/newglst.csv')
#newgeno<-read.csv('c:/Users/jiaxin/Desktop/??????ϵ????/newfrge/snpm.csv')

for(gene in unique(newglst$Gene_Symbol)){
  #for(gene in unique(newglst$Gene_Symbol)){
  #gene<-as.numeric(gene)
  #print(gene)
  gene<-as.character(gene)
  pergene <- fourierExpansion(gene,newgeno,newglst,newmap,rng)
  colnames(pergene)<-paste('gene',gene,seq(1:dim(pergene)[2]),sep = '_')
  geno_expn<-cbind(geno_expn,pergene)
  print(paste('gene',gene,'done!',sep = ' '))
  
}

#write.csv(geno_expn,'c:/Users/jiaxin/Desktop/??????ϵ????/sorteddata/expanded_genotype.csv',
#          row.names = F)









#new defined position
fourierExpansion <-function(gene_idx,geno,gene_list,snp_map,rng){
  #missing(geno)
  snpName<-as.vector( snp_map[,2] )
  snp.position<-as.vector( snp_map[,4] )
  snp.chrom=snp_map[,1]
  
  #gene.chrom<-as.vector( gene_list$Chromosome )
  gene.symbol<-as.vector( gene_list$Gene_Symbol)
  
  #gene.start<-subset(gene_list,Gene_Symbol==gene_idx)$Start
  #gene.end<-subset(gene_list,Gene_Symbol==gene_idx)$End
  #gene.end<-as.vector(gene_list$End
  geno<-as.data.frame(geno)
  
  #sub.idx<-as.character(subset(newposi,Genenames==gene_idx)$SNPnames)
  sub.idx<-glist%in%gene_idx
  sub.idx<-names(geno)[sub.idx]
  #sub.idx<-as.numeric(sub.idx)
  #sub.idx<-snp.chrom == gene_idx
  #sub.idx<-sub.idx & (snp.position>(gene.start -rng) )
  #sub.idx<-sub.idx & (snp.position<(gene.end+rng) )
  #scope<-(snp.position>(gene.start -rng) )&&(snp.position<(gene.end+rng) )
  #if(sum(sub.idx)>0){
  #x<-as.matrix(geno[,sub.idx])
  #  print('get the right locis')
  #}else{sub.idx<-snp.chrom == gene_idx
  #}
  x<-as.matrix(geno[,sub.idx])
  #pos<-as.character(subset(newposi,Genenames==gene_idx)$Position)
  pos<-posivec[sub.idx]
  pos<-as.numeric(pos)
  x<-as.data.frame(x)
  for(i in 1:dim(x)[2]){
    x[,i]<-as.numeric(as.character(x[,i]))
  }
  x<-as.matrix(x)
 
  nsnps <-dim(x)[2]
  
  gil=list();
  # if the number of snp in the gene is over 3, expansion occurs
  if(nsnps>1)
  {
    if ( is.null(pos) )
    {
      pos <-(0:( nsnps-1) )/(nsnps-1)
    }else {
      idx<-order(pos)
      print(pos)
      print(idx)
      print(colnames(x))
      #if(length(order)!=1){
      x<-x[,idx]
      
      pos<-pos[idx]
      pos<-(pos-pos[1])/(pos[nsnps]-pos[1])
    }
    
    # use PCA to determine the number of basis
    eigenval<-prcomp(x)$sd^2
    sum_eigen=sum(eigenval)
    tmp=0
    n_of_basis=0
    for(i in 1:length(eigenval))
    {
      tmp=eigenval[i]+tmp
      n_of_basis=i;
      if(tmp>=0.8*sum_eigen) break
    }
    n_of_basis=floor(n_of_basis/2)*2+1 
    #make n_of_basis_A the odd number
    # end of setting the basis number
    
    frange <-c(pos[1], pos[length(pos)])
    
    fbasis<-create.fourier.basis(frange,nbasis=n_of_basis)
    phi=eval.basis(pos,fbasis);
    gil$zeta <-t(ginv(t(phi)%*%phi)%*%t(phi)%*%t(x))
  }else gil$zeta  <-x
  
  return(gil$zeta)
  
}






###############==================================================================
