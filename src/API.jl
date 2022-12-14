"""
      findGenes(; genome,
                  ref,
                  hasN::Bool = true,
                  mode::String = "return",
                  fileloc::String = "noPath",
                  k::Int64 = 6,
                  windowsize::Int64 = 0,
                  thr::Int64 = 0,
                  buffer::Int64 = 50,
                  resultVec::Vector{FASTA.Record} = FASTA.Record[])

The main API to find all gene matches in a genome from a reference sequence/reference sequence set.
Returns the approximate matches as FASTA records. The descriptions of the record contain information about the match. 
The format of the description appears as so:

   "Identifier | SED = a | Pos = b:c | thr = d | buffer = e"

Where `SED` indicates how similar the match is to the references, and the lower, the more similar in terms of kmer distance.

...
# Arguments
- `genome`: Either a String or a Vector of strings, where each string is the file location of a `.fasta` file containing the genome.
- `ref`: the file location of a fasta file containing the reference sequence(s). The references should be very similar in length. OR a pre-generated reference kmer frequency dictionary.

# Optional Arguments
- `hasN::Bool = true`: Whether both the fasta files from `ref` and the query `genome` contains any undefined nucleotide `N`. If `false` it runs an alternate algorithm that is alot faster than the one accounting for `N` nts. 
- `mode = "return"`: what the function should return the gene matches in. There are three modes, `return`, `print`, `write`. 
Return mode returns a vector of FASTA records. Print mode simply prints the fasta sequences in the REPL. 
Write mode writes the fasta sequences to a fasta file, of which the location has to be given in the `fileloc` argument
- `fileloc::String = "noPath"`: If the mode argument is `"write"`, this should be the location of a fasta file that the results should be written into.
- `k::Int64 = 6`: the kmer length to use for approximate matching. Generally it should probably remain between 5 - 10 and probably does not have a profound impact on how the matches are found. Check out the pre-print for more information
- `windowsize::Int64 = 0`: the size of the returned matches. It should be the average length of all reference sequences and is computed automatically if the argument is left as 0.
- `thr::Int64 = 0`: the squared euclidian distance threshold for sequences. Out of the context of the algorithm, lower values mean matches have to be more similar to the references. If left as 0, it is automatically computed. Once again the pre-print has information on this argument.
- `buffer::Int64 = 0`: the amount of nucleotides left and right of a matched sequence that should be added to the returned fasta sequences
...

It is recommended to leave all other optional arguments as is, especially the buffer = 50 argument. The purpose of the algorithm is to find approximate matches that can then be BLASTed and aligned.

However, playing with the `thr` argument could return more or less matches every time. 
"""
function findGenes(; #FASTQ, RNA and AA compaibility will be added in the future. Also distance metric may be changed in future
   genome::Any, 
   ref::Any,
   hasN::Bool = true,
   mode::String = "return",
   fileloc::String = "noPath",
   k::Int64 = 6, 
   windowsize::Int64 = 0,
   thr::Float64 = 0.0, 
   buffer::Int64 = 50,
   results::Vector{FASTA.Record} = FASTA.Record[])

   #managing undeclared variables
   k < 4 && println("Such a low k value may not yield the most accurate results")
   if windowsize == 0; windowsize = avgRecLen(ref) end

   #gvariables for other GMA
   RV = fill(0.0,5^k)
   ScaleFactor = 1/(2*k) #maybe this could be an argument as well?
   KD = genKmers(k;withN=true) 
   refKFD = nothing 

   #handling KFD 
   if typeof(ref) == String
      refKFD = genRef(k,ref,KD) #generation of kmer frequency dict
   else #assume its dict 
      refKFD = ref
   end
   refKFV = kfv(refKFD,KD) #generation of kmer frequency dict
   
   #threshold generation (memory intensive when benchmarking)
   if thr == 0.0; thr += findRandThr(ref,refKFD,KD, ScaleFactor=ScaleFactor) end #random threshold, based on seed 1112

   #check which verson to use. this can be optimized as quite alot of vars are unessecarily initialized above (will fix in next patch)
   if !hasN
      return gma_Nless_API(genome=genome,ref=ref,thr=thr,mode=mode,fileloc=fileloc,k=k,windowsize=windowsize,buffer=buffer,results=results)
   end #I need to test this though 

   threshold_buffer_tag = " | thr = "*string(round(thr))*" | buffer = "*string(buffer) # Maybe it would be wise to put which record it is

   #genome mining, for string and vector of strings which is the database.
   if typeof(genome) == String #need to incorporate the eucGma here.
      open(FASTA.Reader, genome) do io
         for rec in io
            fill!(RV,0.0)
            gma(k=k, record = rec, refVec = refKFV,
            windowsize = windowsize, kmerDict = KD,
            path=fileloc, thr=thr, buff=buffer,
            curr_kmer_freq = RV, 
            thrbuff=threshold_buffer_tag,
            mode = mode, resultVec = results,
            ScaleFactor = ScaleFactor)
         end
      end
   else
      for str in genome
         open(FASTA.Reader, str) do io
            for rec in io
               fill!(RV,0.0)
               gma(k=k, record = rec, refVec = refKFV,
               windowsize = windowsize, kmerDict = KD,
               path=fileloc, thr=thr, buff=buffer,
               curr_kmer_freq = RV,
               thrbuff=threshold_buffer_tag,
               mode = mode, resultVec = results,
               ScaleFactor = ScaleFactor)
            end
         end
      end
   end
   if mode == "return"; (return results) end
end

export findGenes

#gotta update the doc so that Nless