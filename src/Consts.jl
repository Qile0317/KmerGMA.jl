using BioSequences, FASTX

const SubSeq = LongSubSeq{DNAAlphabet{4}}
const Seq = LongSequence{DNAAlphabet{4}}
const DnaSeq = Union{Seq, SubSeq}
const DnaBits = Dict{BioSequences.DNA, UInt}
const Kfv = Union{Vector{Int64}, Vector{Float64}}

const NUCLEOTIDE_BITS = Dict{BioSequences.DNA, UInt}(
    DNA_A => unsigned(0),
    DNA_C => unsigned(1),
    DNA_G => unsigned(2),
    DNA_T => unsigned(3),
    DNA_N => unsigned(3) # somehow has worked wonders.
) 

"""
    getSeq(seq::FASTA.Record)

Alias to get the dna longsequence of a fasta record
"""
function getSeq(seq::FASTX.FASTA.Record)
     return FASTX.FASTA.sequence(BioSequences.LongSequence{DNAAlphabet{4}}, seq)
end

export getSeq