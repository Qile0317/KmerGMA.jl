"""
      findGenes(; genome,
                  ref::String,
                  mode::String = "return",
                  fileloc::String = "noPath",
                  k::Int64 = 6,
                  windowsize::Int64 = 0,
                  thr::Int64 = 0,
                  buffer::Int64 = 50,
                  BLAST = false,
                  resultVec::Vector{FASTA.Record} = FASTA.Record[])

The main API to find all gene matches in a genome from a reference sequence/reference sequence set.
Returns the approximate matches as FASTA records. The descriptions of the record contain information about the match. 
The format of the description appears as so:

   "Identifier | SED = a | Pos = b:c | thr = d | buffer = e"

Where `SED`(Squared Euclidean Distance) indicates how similar the match is to the references.

...
# Arguments
- `genome`: Either a String or a Vector of strings, where each string is the file location of a `.fasta` file containing the genome.
- `ref::String`: the file location of a fasta file containing the reference sequence(s). The references should be very similar in length

# Optional Arguments
- `mode = "return"`: what the function should return the gene matches in. There are three modes, `return`, `print`, `write`. 
Return mode returns a vector of FASTA records. Print mode simply prints the fasta sequences in the REPL. 
Write mode writes the fasta sequences to a fasta file, of which the location has to be given in the `fileloc` argument
- `fileloc::String = "noPath"`: If the mode argument is `"write"`, this should be the location of a fasta file that the results should be written into.
- `k::Int64 = 6`: the kmer length to use for approximate matching. Generally it should probably remain between 5 - 10 and probably does not have a profound impact on how the matches are found. Check out the pre-print for more information
- `windowsize::Int64 = 0`: the size of the returned matches. It should be the average length of all reference sequences and is computed automatically if the argument is left as 0.
- `thr::Int64 = 0`: the squared euclidian distance threshold for sequences. Out of the context of the algorithm, lower values mean matches have to be more similar to the references. If left as 0, it is automatically computed. Once again the pre-print has information on this argument.
- `buffer::Int64 = 0`: the amount of nucleotides left and right of a matched sequence that should be added to the returned fasta sequences
...

It is recommended ot leave all other optional arguments as is, especially the buffer = 50 argument. The purpose of the algorithm is to find approximate matches that can then be BLASTed and aligned.

Unfinished docs. The BLAST argument is currently useless.
"""
function findGenes(; #FASTQ, RNA and AA compaibility will be added in the future. Also distance metric may be changed in future
   genome, 
   ref::String,
   mode::String = "return",
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
   if thr == 0.0; thr += findthr(ref,refKFD,KD) end
   threshold_buffer_tag = " | thr = "*string(round(thr))*" | buffer = "*string(buffer) # Maybe it would be wise to put which record it is

   #genome mining
   if typeof(genome) == String
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
   else
      for str in genome
         open(FASTA.Reader, str) do io
            for rec in io
               gma(k=k, record = rec, refVec = refKFV,
               windowsize = windowsize, kmerDict = KD,
               path=fileloc, thr=thr, buff=buffer,
               rv=copy(RV), #from benchmarking it seems copying isnt that slow?
               thrbuff=threshold_buffer_tag,
               mode = mode, resultVec = results)
            end
         end
      end
   end
   if mode == "return"; (return results) end
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
