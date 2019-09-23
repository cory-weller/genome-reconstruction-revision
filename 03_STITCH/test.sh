STITCH(
  regionStart = 10000000,
  regionEnd = 10800000,
  buffer = 50000,
  inputBundleBlockSize = 100,
  method = "pseudoHaploid",
  outputdir = outputdir, chr = mouse_chr, posfile = mouse_posfile, genfile = mouse_genfile, bamlist = mouse_bamlist, K = mouse_K, nCores = n_cores, nGen = mouse_nGen)

# ABOVE WORKS


../STITCH.R \
--regionStart=10000000 \
--regionEnd=10800000 \
--buffer=50000 \
--inputBundleBlockSize=100 \
--method=pseudoHaploid \
--outputdir=./ \
--chr=chr19 \
--posfile=pos.txt \
--genfile=gen.txt \
--bamlist=bamlist.txt \
--K=5 \
--nCores=4 \
--nGen=5

# ABOVE WORKS

../STITCH.R \
--regionStart=10000000 \
--regionEnd=10800000 \
--buffer=50000 \
--inputBundleBlockSize=100 \
--method=pseudoHaploid \
--outputdir=./ \
--chr=chr19 \
--posfile=pos.txt \
--bamlist=bamlist.txt \
--K=5 \
--nCores=4 \
--nGen=5

# ABOVE WORKS

../STITCH.R \
--inputBundleBlockSize=100 \
--method=pseudoHaploid \
--outputdir=./ \
--chr=chr19 \
--posfile=pos.txt \
--bamlist=bamlist.txt \
--K=5 \
--nCores=4 \
--nGen=5

# ABOVE WORKS

../STITCH.R \
--method=pseudoHaploid \
--outputdir=./ \
--chr=chr19 \
--posfile=pos.txt \
--bamlist=bamlist.txt \
--K=5 \
--nCores=4 \
--nGen=5

# ABOVE DOESN'T WORK

# Try below with indoor cages;
# add --sampleNames_file=samples.txt;
# shouldn't work
../STITCH.R \
--method=pseudoHaploid \
--sampleNames_file=samples.txt \
--outputdir=./ \
--chr=2L \
--posfile=2L.short.txt \
--bamlist=bamlist.txt \
--K=5 \
--nCores=4 \
--nGen=5

# Above didn't work as expected; add --inputBundleBlockSize=100

# Try below with indoor cages; might work

../STITCH.R \
--inputBundleBlockSize=100 \
--method=pseudoHaploid \
--sampleNames_file=samples.txt \
--outputdir=./ \
--chr=2L \
--posfile=2L.short.txt \
--bamlist=bamlist.txt \
--K=5 \
--nCores=4 \
--nGen=5

# ABOVE WORKED; add 'more real' parameters with K=32 (even though it's likely higher)

nohup ../STITCH.R \
--inputBundleBlockSize=100 \
--method=pseudoHaploid \
--sampleNames_file=samples.txt \
--outputdir=./ \
--chr=2L \
--posfile=2L.txt \
--bamlist=bamlist.txt \
--K=32 \
--nCores=4 \
--nGen=5 &

# ABOVE RAN TO COMPLETION

# Try with K of 128
nohup ../STITCH.R --inputBundleBlockSize=100 --method=pseudoHaploid --sampleNames_file=samples.txt --outputdir=./testKof128 --chr=2L --posfile=2L.txt --bamlist=bamlist.txt --K=128 --nCores=1 --nGen=5 &

