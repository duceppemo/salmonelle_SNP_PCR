#!/bin/bash


############
#          #
#   Help   #
#          #
############


function print_help() {
    echo "\
Usage: bash primer_finder_bbduk.sh [-h] -f[-m]|-q -p <primers.txt> -o <./report/> [-n 2] <sample_file_1> [<sample_file_n>]
Mandatory flags:
    
    -f      Input file is fasta format. Choose \"-f\" or \"-q\".
            Multiple files are accepted and must be separated by a \"space\" character.
            Multiple files will be concatenated.

    -q      Input file is fastq format. Choose \"-f\" or \"-q\".
            Multiple files are accepted and must be separated by a \"space\" character.
            Multiple files will be concatenated.

    -m      Fastq is paired-end (expect two input fastq files).

    -p      Input primer sequence file to test in fasta format.

    -o      Output folder.
Optional flag:
    
    -h      Print this help.

    -n      Mismatches allowed. Default is 1.

    -t      Number of threads. Default is all.
"
exit 1
}


#http://wiki.bash-hackers.org/howto/getopts_tutorial
# Default values
fasta=0
fastq=0
paired=0
primers=""
output=""
cpu=$(nproc)
mismatch=1

options=':fqmp:o:t:n:h'

while getopts "$options" opt; do
    case "$opt" in
        f) fasta=1;;
        q) fastq=1;;
        m) paired=1;;
        p) primers="$OPTARG";;
        o) output="$OPTARG";;
        t) cpu="$OPTARG";;
        n) mismatch="$OPTARG";;
        h) print_help;;
        \?) echo "Invalid option: -"$OPTARG"" >&2; print_help;;
        :) echo "Option -"$OPTARG" requires an argument." >&2; print_help;;
    esac
done

# Get frags and arguments
shift $(($OPTIND - 1))

# Exit if required arguments are missing
if [ -z "$fasta" ] && [ -z "$fastq" ] || [ -z "$primers" ] || [ -z "$output" ]; then 
    echo "All \"-f or -q\", \"-p\" and \"-o\" options are mandatory"
    print_help
fi

# Check if only one of "-f" or "-q" is being used
if [ "$fasta" -eq 1 ] && [ "$fastq" -eq 1 ]; then
    echo "Only one type of input file can be selected: fasta (\"-f\") or fastq (\"-q\")."
    print_help
fi

#check if -m is used with -q
if [ "$paired" -eq 1 ] && [ "$fastq" -ne 1 ]; then
    echo "The paired-end option \"-m\" has to be used the the fastq option \"-q\""
    print_help
fi

# check if two files supplied with the -m option
if [ "$fastq" -eq 1 ] && [ "$paired" -eq 1 ] && [ "$#" -eq 1 ]; then
    echo "Please supply two fastq files for paired-end data."
    print_help
elif [ "$fastq" ] && [ "$paired" ] && [ "$#" -gt 2 ]; then
    echo "$#"
    echo "Only two fastq input files can be used when using the \"-m\" option."
    print_help
fi

# Check if primer file exists
if [ ! -s "$primers" ]; then
    echo "The primer file provided does not exist."
    print_help
fi

# Check if primer file is in fasta format
primer_line1="$(cat "$primers" | head -n 1)"
primer_char1="${primer_line1:0:1}"
if [ "$primer_char1" != ">" ]; then
    echo "Primer file must be in fasta format."
    print_help
fi

# Check if output folder already exists
if [ -d "$output" ]; then
    echo "Output folder already exists. Output files will be overwritten."
else
    mkdir -p "$output"
fi

# Check if input files were supplied
if [ "$#" -eq 0 ]; then
    echo "You must provide at least one input file."
    print_help
fi

# Check format for each input sequence file

# Array to store the input sequence types
seq_types=()

function check_format()
{
    if [ "${1##*.}" == "gz" ]; then #is it gzipped?
        line1="$(zcat "$1" | head -n 1)" #assuming the first line is a sequence header
    else
        line1="$(cat "$1" | head -n 1)"
    fi
    # echo "$line1"

    # The first character of the first line
    char1="${line1:0:1}"
    # echo "$char1"

    if [ $(grep -E "fastq|fq" <<< "$1") ]; then #is input fastq file? Check file extension fisrt.
        if [ "$char1" == "@" ]; then
            seq_types+=("fastq")
        else
            echo "Invalid fastq file: "$1""
            exit 1
        fi
    else #it's a fasta
        if [ "$char1" = ">" ]; then
            seq_types+=("fasta")
        else
            echo "Invalid fasta file: "$1""
            exit 1
        fi
    fi
}

# List all input files
input_seq_files=()
for f in "${@}"; do
    input_seq_files+=("$f")
done
# echo "${input_seq_files[@]}"

# Check input file(s)
for seq in "${input_seq_files[@]}"; do
    # check if input file exists
    if [ -s "$seq" ]; then
        check_format "$seq"
    else
        echo "$seq does not exist."
        exit 1
    fi
done

# check if all sequences the same type
all_types=($(echo "${seq_types[@]}" | tr " " "\n" | sort | uniq | tr "\n" " "))
# echo "${#all_types[@]}";exit

if [ "${#all_types[@]}" -gt 1 ]; then
    echo "Both fasta and fastq files were detected as input."
    print_help
fi

# Check if number of requested cpu does not exced available CPUs
if [ "$cpu" -gt $(nproc) ]; then
    echo "Requested number of CPUs exceds available CPUs."
    print_help
fi

#check is number of missmatch is >=0 and <=3
if [ "$mismatch" -lt 0 ] || [ "$mismatch" -gt 3 ]; then
    echo "Unauthorized mismatch values."
    print_help
fi


######################
#                    #
#     Resources      #
#                    #
######################


#computer performance
mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
memJava="-Xmx"$mem"g"


##########################
#                        #
#   Log + Dependencies   #
#                        #
##########################


#Date
echo -e "$(date)\n" | tee "${output}"/log.txt
echo -e "User: $(whoami)" | tee -a "${output}"/log.txt
echo -e "Processors: "$cpu"" | tee -a "${output}"/log.txt
echo -e "Memory: "$mem"G" | tee -a "${output}"/log.txt

#pipeline version
echo -e "\nprimer_finer_bbduk.sh version 0.1\n" | tee -a "${output}"/log.txt  # $0

#check if depenencies are installed
#if so, log version

# Java
if hash java 2>/dev/null; then 
    java -version 2>&1 1>/dev/null | grep "java version" | tr -d '"' | tee -a "${output}"/log.txt
else
    echo >&2 "java was not found. Aborting." | tee -a "${output}"/log.txt
    # exit 1
fi

# BBDuk
if hash bbduk.sh 2>/dev/null; then 
    bbduk.sh -v 2>&1 1>/dev/null | grep "version" | tee -a "${output}"/log.txt
else
    echo >&2 "bbduk.sh was not found. Aborting." | tee -a "${output}"/log.txt
    # exit 1
fi

#Get the read length
read_length=$(zcat "$1" \
    | head -n 4 \
    | sed -n '2p' \
    | tr -d "\n" \
    | wc -m)

echo ""$read_length"bp reads detected" | tee -a "${output}"/log.txt


###########
#         #
#   I/O   #
#         #
###########


name="$(basename "$1")" #name without path
sampleName="$(cut -d "_" -f 1 <<< "${name%%.*}")" #nameR1 with the ".fastq.gz"


#############################
#                           #
#   Primer Sequence Files   #
#                           #
#############################


declare -a sizes=()
#make list of primer sizes
while IFS= read -r line; do
    l=$(echo $line | tr -d "\n" | tr -d "\r")
    if [ -z "$l" ]; then
        continue
    fi

    char1="${l:0:1}"

    if [ "$char1" = ">" ]; then  # it's the header
        primer_name="$l"
        # echo $primer_name
    else
        #Primer lengths
        primer_length=$(expr length "$l")
        sizes+=("$primer_length")

        #copy line in a primer seqence in file for each length
        echo -e ""$primer_name"\n"$l"" >> "${output}"/"${sampleName}"_primers_"${primer_length}"-mer.fasta
    fi
done < "$primers"

#Remove duplicates and sort array
eval sizes=($(printf "%q\n" "${sizes[@]}" | sort -u))
# echo "${sizes[@]}"; exit

shortest="${sizes[0]}"


########################
#                      #
#   Counting matches   #
#                      #
########################


# #Create a list of all possible strings to search for according to the number of mismatches requested
# takes 20min to generate all the primers for 2 mismatches, plus repeated ones. Switched to python.
# declare -a alphabet=('A' 'T' 'C' 'G')
# for s in "${sizes[@]}"; do
#     while IFS= read -r line; do
#         l=$(echo $line | tr -d "\n" | tr -d "\r")
#         if [ -z "$l" ]; then
#             continue
#         fi

#         char1="${l:0:1}"
#         if [ "$char1" = ">" ]; then  # it's the header
#             primer_name="$l"
#             # echo $primer_name
#         else
#             len=${#l}
#             last_index=$((len-1))
#             snp=${l:last_index:1}

#             declare -a alphabet=('A' 'T' 'C' 'G')
#             counter=0

#             if [ "$mismatch" -eq 1 ]; then  #one mismatch                
#                 for p in $(seq 1 1 $last_index); do
#                     for q in ${alphabet[@]}; do
#                         let counter=counter+1
#                         new=$(sed "s/./$q/$p" <<< $l)
#                         # echo $new
#                         echo -e ""$primer_name"_"${counter}"\n"$new"" >> "${output}"/"${sampleName}"_primers_"${len}"-mer_clipped.fasta
#                     done
#                 done
#             elif [ "$mismatch" -eq 2 ]; then  #two mismatches
#                 for a in $(seq 1 1 "$last_index"); do    # sed index starts at 1
#                     # echo "a: $a"
#                     for b in "${alphabet[@]}"; do
#                         # echo "b: $b"
#                         new=$(sed "s/./"$b"/"$a"" <<< "$l")
#                         if [ "$a" -lt "$last_index" ]; then
#                             for c in $(seq $((a+1)) 1 "$last_index"); do
#                                 # echo "c: $c"
#                                 for d in "${alphabet[@]}"; do
#                                     # echo "d: $d"
#                                     let counter=counter+1
#                                     newer=$(sed "s/./"$d"/"$c"" <<< "$new")
#                                     # echo "$newer"
#                                     echo -e ""$primer_name"_"${counter}"\n"$newer"" >> "${output}"/"${sampleName}"_primers_"${len}"-mer_clipped.fasta
#                                 done
#                             done
#                         fi
#                     done
#                 done
#             fi
#         fi
#     done < "${output}"/"${sampleName}"_primers_"${s}"-mer.fasta
# done

#Create a list of all possible strings to search for according to the number of mismarch requested
if [ "$mismatch" -gt 0 ] && [ "$mismatch" -le 3 ]; then
    python ~/PycharmProjects/mismatch_creator/mismatch_creator.py \
        -m "$mismatch" \
        -i "$primers" \
        -o "${output}"/primers.fasta
else
    cp "$primers" "${output}"/primers.fasta
fi

#First pass with mismatch allowed
if [ "$fastq" -eq 1 ] && [ "$paired" -eq 1 ]; then  # fastq paried-end
    for s in "${sizes[@]}"; do
        bbduk.sh "$memJava" \
            overwrite=t \
            maskmiddle=f \
            rcomp=t \
            k="$s" \
            hdist=0 \
            in1="$1" \
            in2="$2" \
            ref="${output}"/primers.fasta \
            stats="${output}"/"${sampleName}"_"${s}"-mer.counts
    done
else  # fasta or fastq single end
    for s in "${sizes[@]}"; do
        bbduk.sh "$memJava" \
            overwrite=t \
            maskmiddle=f \
            rcomp=t \
            k="$s" \
            hdist="$mismatch" \
            in="$1" \
            ref="${output}"/primers.fasta \
            stats="${output}"/"${sampleName}"_"${s}"-mer.counts
    done
fi

#combine the counts for all the variants of the same primer
# 1035267_REF_10  1   0.00002%
# 1035267_REF_15  1   0.00002%

# get only useful info from bbduk output
for s in "${sizes[@]}"; do
    cat "${output}"/"${sampleName}"_"${s}"-mer.counts \
        | grep -vE "^#" \
        | cut -f 1,2 \
        > "${output}"/"${sampleName}"_"${s}"-mer.counts.txt
    rm "${output}"/"${sampleName}"_"${s}"-mer.counts
done

# Merge all count
cat "${output}"/"${sampleName}"_??-mer.counts.txt \
    | sort -r \
    > "${output}"/"${sampleName}".counts.txt

rm "${output}"/"${sampleName}"_*-mer.counts.txt
rm "${output}"/"${sampleName}"_primers_*-mer.fasta


if [ "$mismatch" -gt 0 ] && [ "$mismatch" -le 3 ]; then
    python ~/PycharmProjects/mismatch_creator/bbduk_count_merger.py \
        -i "${output}"/"${sampleName}".counts.txt \
        -o "${output}"/"${sampleName}".counts.merged.txt
else
    mv "${output}"/"${sampleName}".counts.txt "${output}"/"${sampleName}".counts.merged.txt
fi


##############
#            #
#   Report   #
#            #
##############


[ -e "${output}"/"${sampleName}"_bbduk.tsv ] && rm "${output}"/"${sampleName}"_bbduk.tsv

#make list of primer sizes
while IFS= read -r line; do
    # remove carriage return
    l=$(echo "$line" | tr -d "\n" | tr -d "\r")

    # Skip blank lines
    if [ -z "$l" ]; then
        continue
    fi

    char1="${l:0:1}"

    if [ "$char1" = ">" ]; then  # it's the header
        primer_name="${l:1}"
        # echo $primer_name
    else
        #get the counts from the bbduk count files
        counts=$(cat "${output}"/"${sampleName}".counts.merged.txt | grep -wF "$primer_name" | cut -f 2)
        [ "$counts" ] || counts=0

        printf "%s\t%s\t%s\n" "$primer_name" "$l" "$counts" | tee -a "${output}"/"${sampleName}"_bbduk.tsv
    fi
done < "$primers"


####################
#                  #
#   In silico PCR  #
#                  #
####################


#Get genotype for all positions
haplotype=()

old_IFS=$IFS  # save the field separator           
IFS=$'\n' # new field separator, the end of line           
for line in $(cat "${output}"/"${sampleName}"_bbduk.tsv); do  # skip header    
    if [ $(grep -F "_REF" <<< "$line") ]; then
        ref_position=$(cut -f 1 <<< $line)
        pos=$(cut -d "_" -f 1 <<< "$ref_position")
        ref_sequence=$(cut -f 2 <<< $line)
        ref_allele="${ref_sequence: -1}"
        ref_count=$(cut -f 3 <<< $line)

    else  # elif [ $(grep -F "_ALT") ]; then
        alt_position=$(cut -f 1 <<< $line)
        alt_sequence=$(cut -f 2 <<< $line)
        alt_allele="${alt_sequence: -1}"
        alt_count=$(cut -f 3 <<< $line)

        if [ "$ref_count" -gt "$alt_count" ]; then
            haplotype+=("$ref_allele")
        elif [ "$ref_count" -lt "$alt_count" ]; then
            haplotype+=("$alt_allele")
        elif [ "$ref_count" -eq "$alt_count" ] && [ "$ref_count" -gt 0 ]; then
            echo "Problem at position "$pos": equal counts for both alleles!"  | tee -a "${output}"/log.txt
            haplotype+=(".")
        else
            echo "Position "$pos" not detected in raw reads" | tee -a "${output}"/log.txt
            haplotype+=(".")
        fi
    fi
done          
IFS=$old_IFS # restore default field separator 

# echo "${haplotype[@]}"

# The haplotype
snpSeq="$(echo "${haplotype[@]}" | tr -d " ")"


#store in primer data in hash
# unset clades
# declare -A clades=( \
#     ["clade1a"]="GGCCCCTGGACACCGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTTCACCCCCTAAAT" \
#     ["clade1b"]="GGCCCCTGGACACCGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTCCACCCCCTAAAT" \
#     ["clade1c"]="GGCCTCTGGACACCGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTTCACCCCCTAAAT" \
#     ["clade2a"]="GGCCCCTGGACACCGCCTGGGACGAGTGTCTAGCGCGCAGGGCGCGTTCACCCCCTAAAT" \
#     ["clade2b"]="GGCCTCTGGACACCGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTCCACCCCCTAAAT" \
#     ["clade3"]="GGCCCCTGGACACCGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTTCACCCCCTAAGG" \
#     ["clade4"]="GGCCCCTGGACACCGTCGGGATCAAGTGTCTATTGCGAAGGGTGGGCTCACTCCCTAAGG" \
#     ["clade5"]="GGCCCCTGGACACCGTCGGGATCAAGTGTCTATTGCGAAGGGTGGGCTCACTCCCTCAGG" \
#     ["clade6"]="AGCACCTGGACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCCTCTCCCTAAGG" \
#     ["clade7a"]="AGCATCTGGACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCCTCTCCCTAAGG" \
#     ["clade7b"]="AGCATCTGGACGCCGTCGGGATTAGATGTCCATTGAGAATGGTGGGCCCTCTCCCTAAGG" \
#     ["clade8"]="AGCATCTGGACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCCTCTCTCTAAGG" \
#     ["clade9a"]="AGCATCTGGACGCCGTCGGGATTAGATATCCATTGAGAAGGGTGGGCCCTCTCCCCACGG" \
#     ["clade9b"]="AGCATCTGGACGCCGTCGGGATTAGATATCCATTGAGACGGGTGGGCCCTCTCCCCACGG" \
#     ["clade10"]="AGCATCTGAACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCTTTTCTTTAAGG" \
#     ["clade11"]="AGCATCTGGACGCCGTCGGAATTAGATGTTCATTGAGAATAGTGGGCCCTCTTCCTAAGG" \
# )

# #Short version with only 16 informative SNPs
# unset clades
# declare -A clades=( \
#     ["1_a"]="GGCCGCCGTGAGTCAA" \
#     ["1_b"]="GGCTGCCGTGAGTCAA" \
#     ["2"]="GGCCGCCGAGAGTCAA" \
#     ["3"]="GGCCGCCGTGAGTCAG" \
#     ["4"]="GGCCGCTGTGAGTCAG" \
#     ["5"]="GGCCGCTGTGAGTCCG" \
#     ["6"]="AGCCGCTGTGAGCCAG" \
#     ["7_a"]="AGCTGCTGTGAGCCAG" \
#     ["7_b"]="AGCTGCTGTGATCCAG" \
#     ["8"]="AGCTGCTGTGAGCTAG" \
#     ["9_a"]="AGCTGCTGTAAGCCAG" \
#     ["9_b"]="AGCTGCTGTACGCCAG" \
#     ["10"]="AGCTACTGTGAGCTAG" \
#     ["11"]="AGCTGCTATGATCCAG" \
#     ["12"]="AACCGCTGTGAGCCAG" \
# )

# 2017-09-06
unset clades
declare -A clades=( \
    ["1_a"]="GGCCCCTGGACACCGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTTCACCCCCTAAAT" \
    ["1_b"]="GGCCCCTGGACACCGCCTGGGTCGAGTGTCTGGTGCGCAGGGCGCGTTCACCCCCTAAAT" \
    ["1_c"]="GGCCCCTGGACACTGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTTCACCCCCTAAAT" \
    ["1_d"]="GGCCCCTAGACACCGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTTCACCCCCTAAAT" \
    ["2_a"]="GGCCCCTGGACACCGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTCCACCCCCTAAAT" \
    ["2_b"]="GGCCCCTGGACACCGCCTGGGTCGAGTATCTAGTGCGCAGGGCGCGTCCACCCCCTAAAT" \
    ["3"]="GGCCCCTGGACACCGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTTCACCCTCTAAAT" \
    ["4_a"]="GGCCCCTGGACACCGCCTGGGACGAGTGTCTAGCGCGCAGGGCGCGTTCACCCCCTAAAT" \
    ["4_b"]="GGCCCCTGGACACCGCCTGGGTCGAGTGTCTAGCGCGCAGGGCGCGTTCACCCCCTAAAT" \
    ["4_c"]="GGCCTCTGGACACCGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTCCACCCCCTAAAT" \
    ["5_a"]="GGCCCCTGGACACCGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTTCACCCCCTAAGG" \
    ["5_b"]="GGCCCCTGGACACCGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTTCACCCCCTAAAG" \
    ["5_c"]="AGCCCCTGGACACCGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTTCACCCCCTAAGG" \
    ["5_d"]="GGCCCCTGGACACCGCCTGGGTCGAGTGTCTAGTGCGAAGGGCGCGTTCACCCCCTAAGG" \
    ["6_a"]="GGCCCCTGGACACCGTCTGGGTCGAGTGTCTATTGCGAAGGGCGCGTTCACCCATAAGG" \
    ["6_b"]="GGCCCCTGGACACCGTCTGGGTCGAGTGTCTATTGCGAAGGGCGCGTTCACCCTTAAGG" \
    ["6_c"]="GGCCCCTGGACACCGTCTGGGTCGAGTGTCTATTGCGAAGGGCGCGTTCACCCCTAAGG" \
    ["6_d"]="GGCCCCTGGACACCGTCTGGGTCGAGTGTCTATTGCGAAGGGCGCGTTCACCCGTAAGG" \
    ["7_a"]="GGCCCCTGGACACCGTCTGGGTCGAGTGTCTATTGCGAAGGGCGCGTTCACCCACTAAGG" \
    ["7_b"]="GGCCCCTGGACACCGTCTGGGTCGAGTGTCTATTGCGAAGGGCGCGTTCACCCTCTAAGG" \
    ["7_c"]="GGCCCCTGGACACCGTCTGGGTCGAGTGTCTATTGCGAAGGGCGCGTTCACCCCCTAAGG" \
    ["7_d"]="GGCCCCTGGACACCGTCTGGGTCGAGTGTCTATTGCGAAGGGCGCGTTCACCCGCTAAGG" \
    ["8_a"]="GGCCCCTGGACACCGTCGGGATCAAGTGTCTATTGCGAAGGGCGCGCTCACTCCCTAAGG" \
    ["8_b"]="GGCCCCTGGACACCGTCGGGATCAAGTGTCTATTGCGAAGGGCGGGCTCACTCCCTAAGG" \
    ["9_a"]="GGCCCCTGGACACCGTCGGGATCAAGTGTCTATTGCGAAGGGTGGGCTCACTCCCTAAGG" \
    ["9_b"]="GGCCCCTGGACACCGTCGGGATCAAGTGTTTATTGCGAAGGGTGGGCTCACTCCCTAAGG" \
    ["9_c"]="GGCCCCTGGACACCGTCGGGATCAAGTGTCTATTGAGAAGGGTGGGCTCACTCCCTAAGG" \
    ["9_d"]="GGCCCCTGGACACCGTCGGGATCAAGTGTCTATTGCGAAGGATGGGCTCACTCCCTAAGG" \
    ["9_e"]="GGCCCCTGGACACCGTCGGGATCAAGTGTCTATTGCGAAGGGTGGGCCCACTCCCTAAGG" \
    ["9_f"]="GGCCCCTGGACACCGTCGGGATCAAGTGTCTATTGCGAAGGGTGGGCTCACTCCTTAAGG" \
    ["10"]="GGCCCCTGGACACCGTCGGGATCAAGTGTCTATTGCGAAGGGTGGGCTCACTCCCTCAGG" \
    ["11_a"]="AGCACCTGGACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCCTCTCCCTAAGG" \
    ["11_b"]="AGCACCTGGACGCCGTCGGGATTAAATGTCCATTGAGAAGGGTGGGCCCTCTCCCTAAGG" \
    ["12_a"]="AGCACCTGGACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTAGGCCCTCTCCCTAAGG" \
    ["12_b"]="AGCACCTGGACGCCGTCGGGATTAGATGTCCATTGAGAAGGGCAGGCCCTCTCCCTAAGG" \
    ["13_a"]="AGCATCTGGACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCCTCTCCCTAAGG" \
    ["13_b"]="AGCATCTGGACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCCTCTCCCTAAAG" \
    ["13_c"]="AGCATCTGGACGCCGTCGGGATTAGATATCCATTGAGAAGGGTGGGCCCTCTCCCTAAGG" \
    ["13_d"]="AGCATCTGGACGCCGTCGGGATTGGATGTCCATTGAGAAGGGTGGGCCCTCTCCCTAAGG" \
    ["13_e"]="AGCATCTGGACGCCGTCGGGATTAGATGTCCATTGAGAAGGGCGGGCCCTCTCCCTAAGG" \
    ["14"]="AGCATCTGAACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCTTCTCCCTAAGG" \
    ["15_a"]="AGCATCTGGACGCCGTCGGGATTAGATGTCCATTGAGAATGGTGGGCCCTCTCCCTAAGG" \
    ["15_b"]="AGCATCTGGACGCCGTCGGGATTAGATATCCATTGAGAAGGGTGGGCCCTCTCCCTAAGG" \
    ["16_a"]="AGCATCTGGACGCCGTCGGGATTAGATGTCCATTGAGAATAGTGGGCCCTCTCCCTAAGG" \
    ["16_b"]="AGCATCTGGACGCCGTCGGAATTAGATGTTCATTGAGAATAGTGGGCCCTCTTCCTAAGG" \
    ["17_a"]="AGCATCTGGACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCCTCTCCTTAAGG" \
    ["17_b"]="AGCATCTGGACGCCGTCGGGATTAGATGTCCATTGAGAATGGTGGGCCCTCTCCTTAAGG" \
    ["18"]="AGCATCTGGACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCCTCTCTCTAAGG" \
    ["19"]="AGCATCTGGACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCCTCTCTTTAAGG" \
    ["20_a"]="AGCATCTGGACGCCGTCGGGATTAGATATCCATTGAGAAGGGTGGGCCCTCTCCCCACGG" \
    ["20_b"]="AGCATCTGGACGCCGTCGGGATTAGATATCCATTGAGACGGGTGGGCCCTCTCCCCACGG" \
    ["20_c"]="AGCATCTGGACGCCGTCGGGATTAGATATCCATTGAGACGGGTGGGCCCTCTCTCCACGG" \
    ["21_a"]="AGCATCTGAACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCTTTTCTTTAAGG" \
    ["21_b"]="AGCATCTGAACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCTTTTCTCTAAGG" \
    ["21_c"]="AGCATCTGAACGCCGTCGGGATTAGATGTCCGTTGAGAAGGGTGGGCCTTTTCTTTAAGG" \
    ["21_d"]="AGCATCTGAACGCCGTCGGGATTAGACGTCCATTAAGAAGGGTAGGCCTTTTCTTTAAGG" \
    ["21_e"]="AGCATCTGAACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCTTTTCTCTAAGG" \
    ["21_f"]="AGCATCTGGACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCCTCTCCCTAAGG" \
    ["21_g"]="AGCATCTGAACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGCCTTCTCTCTAAGG" \
    ["21_h"]="AGCATCTGAACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTGGGTCTTTTCTTTAAGG" \
    ["22_a"]="AGCATCTGGACGCCGTCGGAATTAGATGTTCATTGAGAATGGTGGGCCCTCTTCCTAAGG" \
    ["22_b"]="AGCACCTGGACGCCGTCGGGATTAGATGTCCATTGAGAAGGGTAGGCCCTCTCCCTAAGG" \
    ["23"]="AACACCTGGGCGCCGTCGGGATTAGACGTCCGTTAAGAAGGGTAGACCCTCTCCCTAAGG" \
    ["24_a"]="AATACTCAGGTGTTTTAGAGATTAGACGCCCGTTAAAAAGGATAGACCCTCTCCCTAAGG" \
    ["24_b"]="AATACTCAGGCGCTTTAGAGATTAGACGCCCGTTAAAAAGGATAGACCCTCTCCCTAAGG" \
    ["25"]="GACCCCTGGACACCGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTCCACCCCCTAAAT" \
    ["26"]="GGCCTCTGGACACCGCCTGGGTCGAGTGTCTAGTGCGCAGGGCGCGTTCACCCCCTAAAT" \
)

#Identify the clade from the initial array data
for j in "${!clades[@]}"; do
    t=$(grep -E "$snpSeq" <<< "${clades["$j"]}")
    [ "$t" ] && result="$j" && echo "Found clade "$result""
done

#print the clade
if [ -n "$result" ]; then
    echo -e "Isolate "$sampleName" is part of clade "$result""
    touch "${output}"/"$result"
    echo -e "$snpSeq" | tee ""${output}"/"$sampleName"_clade_"${result}".id"
else
    echo -e "Isolate "$sampleName" is not part of any known clade"
    touch "${output}"/unknown_clade
    echo -e "$snpSeq" | tee ""${output}"/"$sampleName"_clade_unknown.id"
fi
