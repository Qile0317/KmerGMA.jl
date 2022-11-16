"""
function findGenes(;genome::FASTX.FASTA.Reader,
   ref::FASTX.FASTA.Reader, fileloc::String,
   k::Int64 = 6, windowsize::Int64 = 0,
   print::Bool = false, thr::Int64 = 0,
   buffer::Int64 = 50, BLAST = true)

The main API to find all matches. Its impportant that the genome to be scanned is a reader object.

The ref argument is a reader of reference FASTA sequences. fileloc is simply the location of a blank file to read the matches to.

- I def need to incorporate option to choose if a plot is wanted instead of SEDs
"""
function findGenes(;genome::FASTX.FASTA.Reader,
   ref::FASTX.FASTA.Reader, fileloc::String,
   k::Int64 = 6, windowsize::Int64 = 0,
   print::Bool = false, thr::Int64 = 0,
   buffer::Int64 = 50, BLAST = true)

   #managing undeclared variables
   k < 4 || error("try a higher value like 6. It is most likely more accurate") #adressing k-value
   windowsize == 0 && (windowsize = avgRecLen(refSeqs,true)) #finding of adequate windowsize
   thr == 0 && (thr = Float64(findthr(ref,refKFV,KD))) #SED threshold estimation

   #variables needed for the GMA
   KD = genKmers(k;withN=true) #generation of kmer dictionary
   refKFV = kfv(genRef(k,ref,KD)) #generation of kmer frequency vector (dict form)
   RV = fill(0.0,5^k)
   threshold_buffer_tag = " | thr = "*string(thr)*" | buffer = "*string(buff)

    #genome mining
    for record in genome
       gma(k=k, record = record, refVec =refKFV,
       windowsize = windowsize, kmerDict = KD,
       path=fileloc, thr=thr, buff=buffer,
       rv=copy(RV), #from benchmarking it seems copying isnt that slow?
       thrbuff=threshold_buffer_tag)
    end
    close(ref)
    close(genome)
    #if BLAST
      #blastseq = eqBLAST(output) #blasting
      #return blastseq #unfinished
   #end
end

export findGenes


#scanning all of vicpac. i should put a progress bar.
#@time writeQueryMatch(6,open(FASTA.Reader,VicPac),V3NRef,genKmers(6,withN=true),350,289,50,"VicPacScan/vicpacscan.fasta")
#@time didnt work but this run started at around 7:06 and ended at 7:23 so it took 17 minutes ish.
