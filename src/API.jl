"""
    findGenes(genome::FASTX.FASTA.Reader,
              ref::FASTX.FASTA.Reader,
              fileloc::String = "";
              k::Int64 = 6,
              windowsize::Int64 = 0,
              print::Bool = false,
              thr::Int64 = 0,
              buffer::Int64 = 50,
              BLAST = true)

The main API to find all matches. Its impportant that the genome to be scanned is a reader object.

The ref argument is a reader of reference FASTA sequences. fileloc is simply the location of a blank file to read the matches to.

- I def need to incorporate option to choose if a plot is wanted instead of SEDs
"""
function findGenes(genome::FASTX.FASTA.Reader,
   ref::FASTX.FASTA.Reader, fileloc::String = "";
   k::Int64 = 6, windowsize::Int64 = 0,
   print::Bool = false, thr::Int64 = 0,
   buffer::Int64 = 50, BLAST = true) #printing is ac pretty dumb

   #managing undeclared variables
   k < 4 || error("try a higher value like 6. It is most likely more accurate") #adressing k-value

   KD = genKmers(k;withN=true) #generation of kmer dictionary

   refKFD = genRef(k,ref,KD) #generation of kmer frequency vector (dict form)

   if windowsize == 0 #finding of adequate windowsize
      windowsize = avgRecLen(refSeqs,true)
   end

    if thr == 0 #SED threshold estimation
        thr = findthr(ref,refKFV,KD)
    end

    #genome mining
    writeQueryMatch(k,genome,refKFD,KD,thr,windowsize,buff,fileloc)

    if BLAST
      blastseq = eqBLAST(output) #blasting
      return blastseq #unfinished
   end
end

export findGenes
