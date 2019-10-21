nohup ~/gatk-4.1.3.0/gatk IndexFeatureFile -F /mnt/internal_1/cory/STITCH_test/indoor_cages/stitch.2L.vcf

counts () {
  ###bamName=/mnt/icy_2/indoorAse/01_mappedData/10.test_run.merge.bam
  bamName=${1}
  echo ${bamName}

  bamStem=$( echo $bamName | rev | cut -d "/" -f1  | rev )

  ~/gatk-4.1.3.0/gatk ASEReadCounter \
  -R /mnt/internal_1/cory/IndoorFlyCages/data_external/all_dmel.fasta \
  -O /mnt/ssd/indoorCages/${bamStem}.alleleCounts.csv \
  -I  ${bamName} \
  -V /mnt/internal_1/cory/STITCH_test/indoor_cages/stitch.2L.vcf \
  --intervals 2L:5000000-6000000 \
  --verbosity ERROR
}
export -f counts

nohup parallel --gnu -j1 -a /mnt/icy_2/indoorAse/01_mappedData/bam.list counts &



### libraries
  library(SNPRelate)
  library(data.table)
  library(stringr)
  library(foreach)
  library(doMC)
  registerDoMC(20)
  library(ggplot2)
  library(ggbeeswarm)

### convert stich VCF to GDS
  #snpgdsVCF2GDS(vcf.fn="/mnt/internal_1/cory/STITCH_test/indoor_cages/stitch.2L.vcf.gz",
  #              out.fn="/mnt/internal_1/cory/STITCH_test/indoor_cages/stitch.32.2L.vcf.gds")


  #snpgdsVCF2GDS(vcf.fn="/mnt/internal_1/cory/STITCH_test/indoor_cages/testKof128/stitch.2L.vcf.gz",
  #              out.fn="/mnt/internal_1/cory/STITCH_test/indoor_cages/stitch.128.2L.vcf.gds")

### open all three genofiles
  geno32 <- snpgdsOpen("/mnt/internal_1/cory/STITCH_test/indoor_cages/stitch.32.2L.vcf.gds")
  geno128 <- snpgdsOpen("/mnt/internal_1/cory/STITCH_test/indoor_cages/stitch.128.2L.vcf.gds")
  genoHS <- snpgdsOpen("/mnt/internal_1/cory/STITCH_test/indoor_cages/filtered.all.vcf.gds")


### get SNP tables
  geno32.dt <- data.table(id=read.gdsn(index.gdsn(geno32, path="snp.id")),
                          chr=read.gdsn(index.gdsn(geno32, path="snp.chromosome")),
                          pos=read.gdsn(index.gdsn(geno32, path="snp.position")),
                          allele=read.gdsn(index.gdsn(geno32, path="snp.allele")),
                          build="geno32")

  geno128.dt <- data.table(id=read.gdsn(index.gdsn(geno128, path="snp.id")),
                          chr=read.gdsn(index.gdsn(geno128, path="snp.chromosome")),
                          pos=read.gdsn(index.gdsn(geno128, path="snp.position")),
                          allele=read.gdsn(index.gdsn(geno128, path="snp.allele")),
                          build="geno128")



  genoHS.dt <- data.table(id=read.gdsn(index.gdsn(genoHS, path="snp.id")),
                          chr=read.gdsn(index.gdsn(genoHS, path="snp.chromosome")),
                          pos=read.gdsn(index.gdsn(genoHS, path="snp.position")),
                          allele=read.gdsn(index.gdsn(genoHS, path="snp.allele")),
                          build="genoHS")


  setkey(geno32.dt, chr, pos)
  setkey(geno128.dt, chr, pos)
  setkey(genoHS.dt, chr, pos)

  dim(geno32.dt)
  dim(geno128.dt)
  dim(genoHS.dt)


### merge
  m <- merge(geno32.dt, geno128.dt, suffixes = c("_32", "_128"))
  m <- merge(m, genoHS.dt, suffixes = c("_foo", "_HS"))

  m <- m[allele_32==allele_128 & allele_32==allele]
  m[,numAllele := str_count(allele, "/")]


### pull in+ ASE calls
  fl <- system("ls /mnt/ssd/indoorCages/*.csv", intern=T)

  obs <- foreach(f=fl, .combine="rbind")%dopar%{
    #f <- fl[1]
    print(f)
    temp <- fread(f)

    num <- first(tstrsplit(last(tstrsplit(f, "/")), "\\."))

    setnames(temp, c("contig", "position"), c("chr", "pos"))

    temp[,geno := -1]
    temp[totalCount>15 & totalCount==refCount, geno:=2]
    temp[totalCount>15 & totalCount==altCount, geno:=0]
    temp[refCount>=2 & altCount>=2, geno:=1]

    temp[,sample.id:=paste("4CB-1-", num, sep="")]

    temp[geno!=-1,c("chr", "pos", "geno", "sample.id"), with=F]
  }


### get genotypes
  ### get SNP ids to use based on informative calls form observed data
    obs[,st:=paste(chr, pos, sep="_")]
    obs.sites <- unique(obs$st)
    obs.sites.dt <- data.table(chr=tstrsplit(obs.sites, "_")[[1]], pos=as.numeric(tstrsplit(obs.sites, "_")[[2]]), key="chr,pos")

    setkey(m, chr, pos)
    setkey(obs.sites.dt, chr, pos)

    snp.ids <- merge(m, obs.sites.dt)$id

  ### pull genotypes

    gt32 <- snpgdsGetGeno(gdsobj=geno32, snp.id=snp.ids, snpfirstdim=F, with.id=T)

    gt32.dt <- data.table(gt32=expand.grid(gt32$genotype)[,1],
                          sample.id=rep(gt32$sample.id, length(gt32$snp.id)),
                          snp.id=rep(gt32$snp.id, each=length(gt32$sample.id)))

    gt128 <- snpgdsGetGeno(gdsobj=geno128, snp.id=snp.ids, snpfirstdim=F, with.id=T)
    gt128.dt <- data.table(gt128=expand.grid(gt128$genotype)[,1],
                          sample.id=rep(gt128$sample.id, length(gt128$snp.id)),
                          snp.id=rep(gt128$snp.id, each=length(gt128$sample.id)))


    gtHS <- snpgdsGetGeno(gdsobj=genoHS, snp.id=snp.ids, snpfirstdim=F, with.id=T)
    gtHS.dt <- data.table(gtHS=expand.grid(gtHS$genotype)[,1],
                          sample.id=rep(gtHS$sample.id, length(gtHS$snp.id)),
                          snp.id=rep(gtHS$snp.id, each=length(gtHS$sample.id)))

    setkey(gt32.dt, sample.id, snp.id)
    setkey(gt128.dt, sample.id, snp.id)
    setkey(gtHS.dt, sample.id, snp.id)

    mgt <- merge(gt32.dt, gt128.dt)
    mgt[,sample.id:=gsub(".bam", "", sample.id)]
    mgt[,sample.id:=tstrsplit(sample.id, "_")[[1]]]
    setkey(mgt, sample.id, snp.id)

    mgt <- merge(mgt, gtHS.dt)

    obs.dt <- merge(m, obs)[,c("chr", "pos", "sample.id", "geno", "id"), with=F]
    setnames(obs.dt, "id", "snp.id")

    setkey(obs.dt, snp.id, sample.id)
    setkey(mgt, snp.id, sample.id)

    mgto <- merge(mgt, obs.dt)



### summaries
    mgto.ag <- mgto[,list(concord=c(mean(gt32==geno, na.rm=T),
                                    mean(gt32==geno, na.rm=T),
                                    mean(gtHS==geno, na.rm=T)),
                        class=c("32", "128", "HS"),
                        n=length(geno)),
                     list(sample.id, geno)]


    ggplot(data=mgto.ag, aes(y=concord, x=class, fill=geno)) + geom_boxplot(position="dodge")

    boxplot(concord~class, mgto.ag)
