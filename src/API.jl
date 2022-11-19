"""
   findGenes(; genome::String,
               ref::String,
               mode = "return",
               fileloc::String = "noPath",
               k::Int64 = 6,
               windowsize::Int64 = 0,
               thr::Int64 = 0,
               buffer::Int64 = 50,
               BLAST = false,
               resultVec::Vector{FASTA.Record} = FASTA.Record[],)

The main API to find all matches. Its impportant that the genome and ref arguments are fasta file locations.

FASTQ, RNA and AA compaibility will be added in the future

The ref argument is a reader of reference FASTA sequences. fileloc is simply the location of a blank file to read the matches to.

There are three modes, return, write, print. Write requires a file location string of a fasta file.

Unfinished docs. The BLAST argument is currently useless.
"""
function findGenes(;
   genome::String, 
   ref::String,
   mode = "return",
   fileloc::String = "noPath",
   k::Int64 = 6, 
   windowsize::Int64 = 0,
   thr::Float64 = 0.0, 
   buffer::Int64 = 50,
   BLAST = false, 
   results::Vector{FASTA.Record} = FASTA.Record[])

   #managing undeclared variables
   k < 4 && error("try a higher value like 6. It is most likely more accurate") #adressing k-value
   if windowsize == 0; windowsize = avgRecLen(ref) end

   #variables needed for the GMA
   KD = genKmers(k;withN=true) #generation of kmer dictionary
   refKFD = genRef(k,ref,KD) #generation of kmer frequency dict
   refKFV = kfv(refKFD,KD) #generation of kmer frequency dict
   RV = fill(0.0,5^k)
   if thr == 0.0
      thr += findthr(ref,refKFD,KD)
   end #SED threshold estimation I tried to do it on 1 line before but it didnt work??
   threshold_buffer_tag = " | thr = "*string(round(thr))*" | buffer = "*string(buffer)

   #genome mining
   open(FASTA.Reader, genome) do io
      for rec in io
         gma(k=k, record = rec, refVec = refKFV,
         windowsize = windowsize, kmerDict = KD,
         path=fileloc, thr=thr, buff=buffer,
         rv=copy(RV), #from benchmarking it seems copying isnt that slow?
         thrbuff=threshold_buffer_tag,
         mode = mode, resultVec = results)
      end
   end
   if mode == "return"; (return results) end
    #if BLAST
      #blastseq = eqBLAST(output) #blasting
      #return blastseq #unfinished
   #end
end
export findGenes

findGenes(genome = "test/Loci.fasta", ref = "test/Alp_V_ref.fasta")

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
