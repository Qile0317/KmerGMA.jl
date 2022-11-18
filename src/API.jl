"""
   findGenes(; genome::FASTX.FASTA.Reader,
               ref::FASTX.FASTA.Reader,
               fileloc::String,
               k::Int64 = 6,
               windowsize::Int64 = 0,
               print::Bool = false,
               thr::Int64 = 0,
               buffer::Int64 = 50,
               BLAST = true)

The main API to find all matches. Its impportant that the genome to be scanned is a reader object.

The ref argument is a reader of reference FASTA sequences. fileloc is simply the location of a blank file to read the matches to.

Unfinished. The print and BLAST arguments are currently useless.
"""
function findGenes(;genome::FASTX.FASTA.Reader,
   ref::FASTX.FASTA.Reader, fileloc::String,
   k::Int64 = 6, windowsize::Int64 = 0,
   print::Bool = false, thr::Int64 = 0,
   buffer::Int64 = 50, BLAST = true)

   #managing undeclared variables
   k < 4 && error("try a higher value like 6. It is most likely more accurate") #adressing k-value
   windowsize == 0 && (windowsize = avgRecLen(ref,true)) #finding of adequate windowsize

   #variables needed for the GMA
   KD = genKmers(k;withN=true) #generation of kmer dictionary
   refKFV = kfv(genRef(k,ref,KD),KD)
   RV = fill(0.0,5^k)
   thr == 0 && (thr = Float64(findthr(ref,refKFV,KD))) #SED threshold estimation
   threshold_buffer_tag = " | thr = "*string(round(thr))*" | buffer = "*string(buffer)

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

"""
#old testing code
open(FASTX.FASTA.Reader, "test/Alp_V_ref.fasta") do reference
    open(FASTX.FASTA.Reader, "test/Alp_V_locus.fasta") do target
        KD = genKmers(6,withN=true)
        k = 6
        refKFD = genRef(k,reference,KD) #generation of kmer frequency dict
        refKFV = kfv(refKFD,KD)
        inp = FASTA.Record[]

        test_gma(k=k,record=first(target),
        refVec = refKFV, windowsize = 289,
        kmerDict =KD,
        path = inp,
        thr = 250.0, buff = 20,
        rv= fill(0.0,5^6),
        thrbuff = " test ")

        print(inp)
    end
end
"""

"""
   testFindGenes(;genome::String,
                  ref::String,
                  k::Int64 = 6,
                  windowsize::Int64 = 0,
                  thr::Int64 = 0,
                  buffer::Int64 = 50)

testing function only. Returns Vector of FASTA records.
"""
function testFindGenes(;
   genome::String,
   ref::String,
   k::Int64 = 6,
   windowsize::Int64 = 0, thr::Float64 = 0.0, buffer::Int64 = 50)

   k < 4 && error("try a higher value like 6. It is most likely more accurate") #adressing k-value
   windowsize == 0 && (windowsize = avgRecLen(ref)) #finding of adequate windowsize

   #variables needed for the GMA
   KD = genKmers(k;withN=true) #generation of kmer dictionary
   refKFD = genRef(k,ref,KD) #generation of kmer frequency dict
   refKFV = kfv(refKFD,KD) #generation of kmer frequency dict
   RV = fill(0.0,5^k)
   results = FASTX.FASTA.Record[]
   if thr == 0.0
      thr += findthr(ref,refKFD,KD)
   end #SED threshold estimation I tried to do it on 1 line before but it didnt work??
   threshold_buffer_tag = " | thr = "*string(round(thr))*" | buffer = "*string(buffer)

   open(FASTX.FASTA.Reader, genome) do reader
      for rec in reader
         test_gma(k=k, record = rec, refVec =refKFV,
         windowsize = windowsize, kmerDict = KD,
         path=results, thr=thr, buff=buffer,
         rv=copy(RV), #from benchmarking it seems copying isnt that slow?
         thrbuff=threshold_buffer_tag)
      end
   end
   return results
end
export testFindGenes
