using BioSequences, FASTX

export SubSeq, Seq, DnaSeq, DnaBits, Kfv

const SubSeq = BioSequences.LongSubSeq{DNAAlphabet{4}}
const Seq = BioSequences.LongSequence{DNAAlphabet{4}}
const DnaSeq = Union{Seq, SubSeq}
const DnaBits = Dict{BioSequences.DNA, UInt}
const Kfv = Union{Vector{Int64}, Vector{Float64}}

const JuliaPalette = Dict{String, String}(
    "purple" => "#9358A4",
    "red" => "#CB392E",
    "green" => "#369844",
    "blue" => "#4C64B0"
)

const NUCLEOTIDE_BITS = Dict{BioSequences.DNA, UInt}(
    DNA_A => unsigned(0),
    DNA_C => unsigned(1),
    DNA_G => unsigned(2),
    DNA_T => unsigned(3),
    DNA_N => unsigned(3) # somehow has worked wonders.
) 

export NUCLEOTIDE_BITS

"""
    getSeq(sequence::FASTA.Record)

Alias to get the dna longsequence of a fasta record
"""
function getSeq(sequence::FASTX.FASTA.Record)
     return FASTX.FASTA.sequence(Seq, sequence)
end

export getSeq

getK(KFV::Kfv) = Int(log(4, length(KFV)))
export getK