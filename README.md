## Beatrice Parodi, 7.06.2017
#### Population genetics of killer whales

Before starting the analysis, create a folder with data and results and make sure that all files are in the rigth format

1. Make different lists of bamfiles (all and each population)
Create a path for data and then lists

`DATA=~/Desktop/killerwhales/data`

`ls $DATA/typeC_*.sorted.bam > typeC.bamlist`
`ls $DATA/B1_*.sorted.bam > typeB1.bamlist`
`ls $DATA/B2_*.sorted.bam > typeB2.bamlist`
`ls $DATA/transient_*.sorted.bam > transient.bamlist`
`ls $DATA/resident_*.sorted.bam > resident.bamlist`

2. Create a list with all bamfiles

`ls $DATA/*.sorted.bam > ALL.bamlist`

3. Set the paths

`ANGSD=~/Desktop/install/angsd/angsd`
`NGSTOOLS=~/Desktop/install/ngsTools`
`NGSADMIX=~/Desktop/install/ngsadmix/NGSadmix `
`NGSLD=~/Desktop/install/ngsLD`
 
## Analyse dataset (Qscores & depth)

1. Create a path for reference and print the quality scores and depth (global and per sample)

`REF=~/Desktop/killerwhales/data/unplaced.scaf.fna.gz`

`$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out results/ALL.qc -minMapQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 500 &> /dev/null`

2. Look at the results with

`less -S results/ALL.qc*`

3. Plot to visualise the results with Rscript

`Rscript $NGSTOOLS/Scripts/plotQC.R results/ALL.qc`

`less -S results/ALL.qc.info`

`evince results/ALL.qc.pdf`

4. calculate average depth 

`$ANGSD/angsd -bam ALL.bamlist -doDepth 1 -out results/ALL.depth.Angsd -doCounts 1 -r KB316843.1: -nInd 50`

5. Print individual depth per sample

`$ANGSD/angsd -out Individialdepth -r KB316843.1: -doCounts 1 -dumpCounts 2 -bam ALL.bamlist`

6. R heatmap 

`Individualdepth <-read.table ("results/ALL.depth.Angsd.depthSample")`

`aa=apply(Individualdepth, 1, function(x) x/sum(x))`

`aa=t(aa)`

`plot(aa[1,])`

`bb=aa[,1:15]`

`my_palette <- colorpanel (75, "slategray1","skyblue2", "royalblue4")` 

`heatmap.2(bb,dendrogram="none",trace="none",Rowv=NA,Colv=NA,density.info="none",margins =c(4,4),col= my_palette ,xlab="Depth",ylab="Individual",main="Depth per individual")`

`pdf("./results/Depthmap.pdf")`

`heatmap.2(bb,dendrogram="none",trace="none",Rowv=NA,Colv=NA,density.info="none",margins =c(4,4),col= my_palette ,xlab="Depth",ylab="Individual",main="Depth per individual")`

`dev.off ()`

7. Weighted mean (values per sample that we are going to use for the Gibraltar simulation)

`R`

read as matrix 

`X = as.matrix("results/ALL.depth.Angsd.depthSample", nrow=50, ncol=101)`

`D = 1:101`

calculate weighted depths

`meanDepth = apply(X, 1, function(t) sum(t*D)/sum(t))`

`wMean = round(meanDepth)`

now we only keep values from 20 to 40 (transient and resident populations) and create a file with the depths

`cat(wMean[20:40], file="sampleDepth.txt")`



## Principal component analysis (PCA)

Performing a PCA

1. First, set thresholds:

* setMinDepth 10
* setMaxDepth 200
* r  KB316842.1:
* SNP_pval 1e-3
* minMapQ 20
* minQ 20
* minInd 25 (half of the total)
* -doGeno 32 for posterior probabilities
* -GL 1
* -doCounts 1

2. Run the command to calculate genotype probabilities with

`$ANGSD/angsd -b killerwhales/ALL.bamlist -ref data/unplaced.scaf.fna.gz -out results/ALL -r KB316843.1: -setMinDepth 10 -setMaxDepth 200 -minMapQ 20 -minQ 20 -minInd 25 -SNP_pval 1e-3 -doMajorMinor 1 -doMaf 1 -doGeno 32 -doCounts 1 -GL 1 -doPost 1 &> /dev/null`

3. Unzip our results in binary format 

`gunzip results/ALL.geno.gz`

Have a look at the output with

`less -S results/ALL.mafs.gz`

Look for how many sites we have with 

`NSITES=`zcat results/ALL.mafs.gz | tail -n+2 | wc -l` `

4. Run NGSCovar for covariance matrix

`$NGSTOOLS/ngsPopGen/ngsCovar -probfile results/ALL.geno -outfile results/ALL.covar -nind 50 -call 0 -norm 0 &> /dev/null`

5. Create a cluster and plot in R

`Rscript -e 'write.table(cbind(seq(1,50),rep(1,50),c(rep("B1",8),rep("B2",11),rep("resident",11),rep("transient",10),rep("C",10))), row.names=F, sep=" ", col.names=c("FID","IID","CLUSTER"), file="results/ALL.clst", quote=F)'`

`Rscript $NGSTOOLS/Scripts/plotPCA.R -i results/ALL.covar -c 1-2-3-4 -a results/ALL.clst -o results/ALL.pca.pdf`

`evince results/ALL.pca.pdf`

6. PCA for Antarctic ecotypes and plot

`$ANGSD/angsd -b antarctic.bamlist -ref data/unplaced.scaf.fna.gz -out results/Antarctic -doMajorMinor 1 -doMaf 1 -doGeno 32 -doCounts 1 -GL 1 -doPost 1`

`gunzip results/Antarctic.geno.gz`

`less -S results/Antarctic.mafs.gz`

`NSITES=`zcat results/Antarctic.mafs.gz | tail -n+2 | wc -l``

`$NGSTOOLS/ngsPopGen/ngsCovar -probfile results/Antarctic.geno -outfile results/Antarctic.covar -nind 50 -call 0 -norm 0`

`Rscript -e 'write.table(cbind(seq(1,29),rep(1,29),c(rep("B1",8),rep("B2",11),rep("C",10))), row.names=F, sep=" ", col.names=c("FID","IID","CLUSTER"), file="results/Antarctic.clst", quote=F)'`

`Rscript $NGSTOOLS/Scripts/plotPCA.R -i results/Antarctic.covar -c 1-2-3-4 -a results/Antarctic.clst -o results/Antarctic.pca.pdf`

`evince results/Antarctic.pca.pdf`


##NGSadmix : calculate admixture proportions

1. Create BEAGLE files with ANGSD

`$ANGSD/angsd -b ALL.bamlist -ref data/unplaced.scaf.fna.gz -out results/admix -r KB316843.1: -doMajorMinor 4 -setMinDepth 10 -setMaxDepth 200 -minMapQ 20 -minQ 20 -minInd 25 -doMaf 1 -doCounts 1 -doGlf 2 -GL 1`

2. Set how many ancestral components do we want to test

`K=5`
 
`$NGSADMIX -likes results/admix.beagle.gz -K $K -minMaf 0.01 -minInd 25 -outfiles results/admixall`

3. Plot the results in R

`R`

`r=read.table("results/admixall.qopt")`

`barplot(t(as.matrix(r,nrow=50)),beside=F , col=rainbow(5), xlab= "Individual #" , ylab= "Ancestry")`

And save the plot with 

`pdf("./results/barplot.pdf")`

`barplot(t(as.matrix(r,nrow=50)),beside=F, col=rainbow(5), xlab= "Individual #" , ylab= "Ancestry")`

`dev.off()`

`q()`

## Site frequency spectrum

1. Perform analysis for each population with a loop

`for POP in typeB1 typeB2 resident transient typeC`

`do`

      `echo $POP`

      `$ANGSD/angsd -b $POP.bamlist -ref data/unplaced.scaf.fna.gz -anc data/unplaced.scaf.fna.gz -out results/$POP -setMinDepth 10 -setMaxDepth 200 -GL 1 -doSaf 1 -doCounts 1 `

`done`
 
Have a looking at the output (for instance look at the transient saf.idx file)

`$ANGSD/misc/realSFS print results/transient.saf.idx | less -S`

2. SFS with realSFS

`for POP in typeB1 typeB2 resident transient typeC`

`do`

     `echo $POP`
     `$ANGSD/misc/realSFS results/$POP.saf.idx 2> /dev/null > results/$POP.sfs`

`done`

And have a look at the output (for instance B1 pod)

`cat results/typeB1.sfs`

3. Plot results in R

`Rscript $NGSTOOLS/Scripts/plotSFS.R results/typeB1.sfs results/typeB2.sfs results/typeC.sfs`

`evince results/B1-B2-C.pdf`

4. Generate bootstrapping replicates

`for POP in typeB1 typeB2 resident transient typeC`

`do`

     `echo $POP`
     `$ANGSD/misc/realSFS results/$POP.saf.idx -bootstrap 10  2> /dev/null > results/$POP.boots.sfs`

`done`

Look at results, for instance at B1 pod, with

`cat results/typeB1.boots.sfs`

5. Perform a 2D SFS. Here how to compare all populations with typeC pod

`for POP in typeB1 typeB2`

`do`

     `echo $POP`
     `$ANGSD/misc/realSFS results/$POP.saf.idx results/typeC.saf.idx 2> /dev/null > results/$POP.typeC.sfs`

`done`

Compare typeB1 and typeB2

`$ANGSD/misc/realSFS results/typeB2.saf.idx results/typeB1.saf.idx 2> /dev/null > results/typeB2.typeB1.sfs`

And plot the results in R. For instance

`Rscript $NGSTOOLS/Scripts/plot2DSFS.R results/typeB1.typeC.sfs 8 10`
`Rscript $NGSTOOLS/Scripts/plot2DSFS.R results/typeB2.typeC.sfs 11 10`
`Rscript $NGSTOOLS/Scripts/plot2DSFS.R results/typeB2.typeB1.sfs 11 8` 

6. Perform 3D SFS

`$ANGSD/misc/realSFS results/typeB1.saf.idx results/typeB2.saf.idx results/typeC.saf.idx 2> /dev/null > results/antarctic.sfs`

## FST & PBS

1. per site FST

`$ANGSD/misc/realSFS fst index results/typeB2.saf.idx results/typeB1.saf.idx results/typeC.saf.idx -sfs results/typeB2.typeB1.sfs -sfs results/typeB2.typeC.sfs -sfs results/typeB1.typeC.sfs -fstout results/FSTtypeC.pbs &> /dev/null`

2. look at the output with

`$ANGSD/misc/realSFS fst print results/FSTtypeC.pbs.fst.idx | less -S`

3. Compute a sliding window analysis

`$ANGSD/misc/realSFS fst stats2 results/FSTtypeC.pbs.fst.idx -win 1 -step 1 > results/FSTtypeC.pbs.txt`

4. look at the output with

`less -S results/FSTpop.pbs.txt`

5. and then calculate the others

` $ANGSD/misc/realSFS fst index results/typeB2.saf.idx results/typeB1.saf.idx -sfs results/typeB2.typeB1.sfs -fstout results/typeB2.typeB1`

`$ANGSD/misc/realSFS fst stats2 results/typeB2.typeB1.fst.idx -win 1 -step 1 > results/typeB2.typeB1.fst.txt`

`less -S results/typeB2.typeB1.fst.txt`

6. PBS

`$ANGSD/misc/realSFS fst index results/typeB2.saf.idx results/typeB1.saf.idx -sfs results/typeB2.typeB1.sfs -fstout results/typeB2.typeB1`

## Nucleotide diversity

1. Calculate allele frequency posterior probabilities and statistic associates 

`for POP in typeB1 typeB2 typeC`

`do`
    
      `echo $POP`
      `$ANGSD/angsd -b $POP.bamlist -ref data/unplaced.scaf.fna.gz -anc data/unplaced.scaf.fna.gz -out results/nd$POP -nInd 29-doCounts 1 -GL 1 -doSaf 1 -doThetas 1 -pest results/$POP.sfs`

`done`

2. Index these files

`for POP in typeB1 typeB2 typeC`

`do`

    `echo $POP`
     `$ANGSD/misc/thetaStat print results/nd$POP.thetas.gz`
     `$ANGSD/misc/thetaStat do_stat results/nd$POP.thetas.gz -win 1 -step 1 - outnames results/$POP.thetas`

`done`

3. And look at the output (for instance of B2) with

`less -S results/typeB2.thetas.pestPG`
 	

## ngsLD, linkage disequilibrium

`NGSLD=/home/beatrice/Desktop/install/ngsLD`

1. Reference file with the site coordinates with

#$SIM_DATA/testA.geno | perl -s -p -e 's/0 0/0/g; s/(\w) \1/2/g; s/\w \w/1/g; $n=s/2/2/g; tr/02/20/ if($n>$n_ind/2)' -- -n_ind=$N_IND | awk '{print "chrSIM\t"NR"\t"$0}' | gzip -cfn --best > testLD_T.geno.gz

#zcat testLD_T.geno.gz | awk 'BEGIN{cnt=1} pos > 10000 {pos=0; cnt++} {pos+=int(rand()*1000+1); print $1"_"cnt"\t"pos}' > testLD.pos

2. ANGSD with -doGeno 32 

 `for POP in B1 B2 C resident transient`

    `do`
    `echo $POP`

    `$ANGSD/angsd -b ALL.bamlist -ref $REF -out results/LD$POP -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 1 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-4 -doGeno 32 -doPost 1`

`done`

`zcat results/LD$POP.geno.gz > results/LD$POP.geno`

`NSITES=$'((zcat results/LD$POP.mafs.gz | wc -l))'`

`zcat results/LD$POP.mafs.gz | cut -f 1,2 | tail -n +2 > results/$POP.pos.txt`

`$NGSLD/ngsLD --verbose 1 --n_ind 50 --n_sites $NSITES --geno results/$POP.geno --probs --pos results/$POP.pos.txt --rnd_sample 0.05 --max_kb_dist 1000 > results/$POP.ld`

`python Fit_Exp.py --input_type FILE --input_name results/$POP.ld --data_type r2GLS --plot`


## D statistics, abbababa test

1. symbolic links to angsd

`ln -s ../install/angsd/angsd/angsd ANGSD`
`ln -s ../install/angsd/angsd/R/estAvgError.R DSTAT`


2. concatenate lists of bam files

`cat typeB1.bamlist typeB2.bamlist typeC.bamlist resident.bamlist > bamResident.filelist`
`cat typeB1.bamlist typeB2.bamlist typeC.bamlist transient.bamlist > bamTransient.filelist`
`cat typeB1.bamlist typeC.bamlist transient.bamlist resident.bamlist > bamTestTrans.filelist`
`cat typeB2.bamlist typeC.bamlist transient.bamlist resident.bamlist > bamTestTrans2.filelist`

3. run Dstat counts of ABBA,BABA alleles combinations on (((B1 B2) C) Resident) tree


`nohup ./ANGSD -doAbbababa2 1 -doCounts 1 -sizeFile resident.size -bam bamResident.filelist -out bamResident.Filter.Angsd -useLast 1 -blockSize 10000 -p 3 -minQ 20 -minmapQ 20 -r KB316843.1: > residentFilter.print &`

#look at output 
#watch tail residentFilter.print 
#calculate the Dstat in R

`Rscript DSTAT angsdFile=bamResident.Filter.Angsd nameFile=resident.name out=resident.R`

4. run Dstat counts of ABBA,BABA alleles combinations on (((B1 B2) C) Transient) tree

`nohup ./ANGSD -doAbbababa2 1 -doCounts 1 -sizeFile transient.size -bam bamTransient.filelist -out bamTransient.Filter.Angsd -useLast 1 -blockSize 10000 -p 3 -minQ 20 -minmapQ 20 -r KB316843.1: > transientFilter.print &`

`Rscript DSTAT angsdFile=bamTransient.Filter.Angsd nameFile=transient.name out=transient.R`

##5. run Dstat counts of ABBA,BABA alleles combinations on (((B1 C) Transient) Resident) tree

##`nohup ./ANGSD -doAbbababa2 1 -doCounts 1 -sizeFile testTrans.size -bam bamTestTrans.filelist -out bamTestTrans.Filter.Angsd -useLast 1 -blockSize 10000 -p 3 -minQ 20 -minmapQ 20 -r KB316843.1: > testTransFilter.print &`

##`Rscript DSTAT angsdFile=bamTestTrans.Filter.Angsd nameFile=testTrans.name out=testTrans.R`

##6. run Dstat counts of ABBA,BABA alleles combinations on (((B1 C) Transient) Resident) tree

##`nohup ./ANGSD -doAbbababa2 1 -doCounts 1 -sizeFile testTrans2.size -bam bamTestTrans2.filelist -out bamTestTrans2.Filter.Angsd -useLast 1 -blockSize 10000 -p 3 -minQ 20 -minmapQ 20 -r KB316843.1: > testTransFilter2.print &`

##`Rscript DSTAT angsdFile=bamTestTrans2.Filter.Angsd nameFile=testTrans2.name out=testTrans2.R`

7. ancient fasta ((B1 B2) C) ancestral) 

`cat typeB1.bamlist typeB2.bamlist typeC.bamlist > bamAncient.filelist`

`ANC=~/Desktop/killerwhales/data/ancestral.fa` 

`nohup ./ANGSD -doAbbababa2 1 -doCounts 1 -sizeFile ~/Desktop/Dstat/ancestral.size  -bam bamAncient.filelist -anc $ANC -out ancestral.Angsd -useLast 0 -blockSize 10000 -p 3 -minQ 20 -minmapQ 20 -r KB316843.1: > ancestral.print &`

8. Look at the output with `watch tail ancestral.print`

9. Calculate the Dstat in R

`Rscript DSTAT angsdFile=ancestral.Angsd nameFile=ancestral.name out=ancestral.R`

# MSc-EEC_Population-genomics-of-killer-whales_2017
Pipeline for MSc project (EEC), 2017
