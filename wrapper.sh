export cpu=$(nproc)  # All CPUs
export mem="$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100))" #85% of total available memory

export data="/media/6tb_raid10/data/salmonella_wallid/fastq/Enteritidis"
export output=""${HOME}"/Desktop/salmonella_test"

export prog=""${HOME}"/prog/genomescope"

function samID ()
{
    r1="$1"
    r2=$(sed 's/_R1_/_R2_/' <<< "$r1")

    bash ~/scripts/salmonella_SNP_PCR.sh \
        -t 4 \
        -q -m \
        -p '/media/6tb_raid10/data/salmonella_clades/SE_clades_primers_21-mer.fasta' \
        -o "$output" \
        "$r1" "$r2"
}

export -f samID

find "$data" -type f -name "*_R1_*" |
parallel    --bar \
            --env samID \
            -j 12 \
            'samID {}'
        


function genomeStats()
{
    r1="$1"
    r2=$(sed 's/_R1_/_R2_/' <<< "$r1")

    name="$(basename "$1")" #name without path
    sampleName="$(cut -d "_" -f 1 <<< "${name%%.*}")" #nameR1 with the ".fastq.gz"

    zcat "$r1" "$r2"\
        | jellyfish count \
            -s "$mem" \
            -m 21 \
            -C \
            -t 8 \
            -o "${output}"/"${sampleName}".jf \
            -L 3 \
            /dev/stdin

    jellyfish histo \
        -t 8 \
        -o "${output}"/"${sampleName}".histo \
        "${output}"/"${sampleName}".jf

    read_length=$(zcat "$1" \
        | head -n 4 \
        | sed -n '2p' \
        | tr -d "\n" \
        | wc -m)

    Rscript "${prog}"/genomescope.R \
        "${output}"/"${sampleName}".histo \
        21 \
        "$read_length" \
        "${output}"/"${sampleName}"
}

export -f genomeStats

find "$data" -type f -name "*_R1_*" |
parallel    --bar \
            --env genomeStats \
            --env output \
            --env mem \
            -j 6 \
            'genomeStats {}'

for i in $(find /home/bioinfo/Desktop/salmonella_test -type f -name "model.txt"); do
    name="$(dirname "$i")"
    sampleName=$(basename "$name")

    coverage=$(cat "$i" \
                | grep -E "^kmercov" \
                | tr -s " " \
                | cut -d " " -f 2)
    cov=$(printf '%.0f' "$coverage")
    echo -e ""$sampleName"\t"$cov"" >> coverages.txt
done

