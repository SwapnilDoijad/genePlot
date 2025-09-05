#!/bin/bash

# Initialize variables
    data_dir=""
    output_dir=""
    cpus=1  # Default value for CPUs

# Function to display usage
    usage() {
        echo "Usage: $0 -d <data_directory> -o <output_directory> [-c <cpus>]"
        echo "  -d, --data_dir   <directory>     Input data directory"
        echo "  -o, --output_dir  <directory>     Output directory"
        echo "  -c, --cpus        <number>       Number of CPUs to use (optional, default: 4)"
        exit 1
    }

# Parse command-line arguments
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            -d|--data_dir)
                data_dir="$2"
                shift 2
                ;;
            -o|--output_dir)
                output_dir="$2"
                shift 2
                ;;
            -c|--cpus)
                cpus="$2"
                shift 2
                ;;
            -h|--help)
                usage
                ;;
            *)
                echo "Unknown option: $1"
                usage
                ;;
        esac
    done

# Validate inputs
    if [[ -z "$data_dir" || -z "$output_dir" ]]; then
        echo "Error: Both data directory and output directory are required."
        usage
    fi

# Debugging: Print the inputs
    echo "---------------------------------------------------------------------"
    echo "Data Directory: $data_dir"
    echo "Output Directory: $output_dir"
    echo "CPUs: $cpus"
    echo "---------------------------------------------------------------------"

## step-0: preparations
    sx=$(ls $data_dir | head -1 | sed -n 's/.*\.//p')
    list=$(find $data_dir -maxdepth 1 -type f -name "*.${sx}" -exec basename {} .${sx} \;)
    log() {
        echo "$(date +"%Y-%m-%d %H:%M:%S"): $1"
    }

## step-1: gff2fasta
    mkdir -p ${output_dir}/fasta > /dev/null 2>&1
    for i in $list; do
        if [ ! -f ${output_dir}/fasta/${i}.fasta ] ; then
            log "Processing GFF file: $i"

            python scripts/gff2fasta.py \
            --gff ${data_dir}/${i}.${sx} \
            --output ${output_dir}/fasta/${i}.fasta
        fi
    done

## step-2: extract CDS
    mkdir -p ${output_dir}/cds > /dev/null 2>&1
    for i in $list; do
        if [ ! -f ${output_dir}/cds/${i}.cds.fasta ] ; then
            log "Extracting CDS for: $i"

            python scripts/extract_cds_for_blast.py \
            -i ${data_dir}/${i}.${sx} \
            -o ${output_dir}/cds/${i}.cds.fasta
        fi
    done

## step-3: run mash sketch
    mkdir -p ${output_dir}/mash > /dev/null 2>&1
    for i in $list; do
        if [ ! -f ${output_dir}/mash/$i.mash_sketch.msh ] ; then
            log "Running mash sketch for: $i"

            ( mash sketch \
            -s 1000 -p $cpus \
            -o ${output_dir}/mash/$i.mash_sketch.msh \
            ${output_dir}/fasta/$i.fasta ) > /dev/null 2>&1
        fi
	done

## step-4: run mash distance
    mkdir -p ${output_dir}/distance > /dev/null 2>&1
    for i in $list; do
	for i2 in $list; do 
        if [ ! -f ${output_dir}/distance/$i.$i2.dist.tsv ] ; then
            log "Running mash distance for: $i"

            mash dist \
            -p $cpus \
			${output_dir}/mash/$i.mash_sketch.msh \
			${output_dir}/mash/$i2.mash_sketch.msh \
			> ${output_dir}/distance/$i.$i2.dist.tsv
        fi
	done
	done

## step-5: calculate shared genome content
	mkdir -p ${output_dir}/SGC_ANI/tmp
	for i in $list ; do
	for i2 in $list ; do
	if [ ! -f ${output_dir}/SGC_ANI/tmp/$i.$i2.txt ] ; then
		log "RUNNING: calculating SGC_ANI $i"

		awk '{print $5}' ${output_dir}/distance/$i.$i2.dist.tsv \
		| awk -F'/' '{print $1*100/1000 }' | tr ' ' '\n' \
		> ${output_dir}/SGC_ANI/tmp/$i.$i2.SGC.txt

		awk '{print (1-$3)*100}' ${output_dir}/distance/$i.$i2.dist.tsv \
		| tr ' ' '\n' > ${output_dir}/SGC_ANI/tmp/$i.$i2.ANI.txt

		paste ${output_dir}/SGC_ANI/tmp/$i.$i2.SGC.txt \
		${output_dir}/SGC_ANI/tmp/$i.$i2.ANI.txt \
		> ${output_dir}/SGC_ANI/tmp/$i.$i2.SGC_ANI.txt

		awk '{print $1*$2/100}' ${output_dir}/SGC_ANI/tmp/$i.$i2.SGC_ANI.txt \
		> ${output_dir}/SGC_ANI/tmp/$i.$i2.txt

		cat ${output_dir}/SGC_ANI/tmp/$i.$i2.txt \
		>> ${output_dir}/SGC_ANI/tmp/$i2.txt
	fi
	done
    done

## step-6: summarise SGC_ANI results
	if [ ! -f ${output_dir}/SGC_ANI/SGC_ANI.tab ]; then 
		printf "%s\n" $list > ${output_dir}/SGC_ANI/SGC_ANI.txt
        for i in $list; do 
			log "SUMMARISING: mash SGC for $i"

            paste \
            ${output_dir}/SGC_ANI/SGC_ANI.txt \
            ${output_dir}/SGC_ANI/tmp/$i.txt \
            > ${output_dir}/SGC_ANI/SGC_ANI.tmp

            mv ${output_dir}/SGC_ANI/SGC_ANI.tmp ${output_dir}/SGC_ANI/SGC_ANI.txt
        done

        printf "%s\t" $list > ${output_dir}/SGC_ANI/tmp/header.tmp
        sed -i '1s/^/\t/' ${output_dir}/SGC_ANI/tmp/header.tmp
        sed -i -e '$a\' ${output_dir}/SGC_ANI/tmp/header.tmp

        cat ${output_dir}/SGC_ANI/tmp/header.tmp ${output_dir}/SGC_ANI/SGC_ANI.txt \
		> ${output_dir}/all.SGC_ANI.tmp
        
		sed -i 's/ /\t/g' ${output_dir}/all.SGC_ANI.tmp

        mv ${output_dir}/all.SGC_ANI.tmp ${output_dir}/SGC_ANI/SGC_ANI.tab
        rm ${output_dir}/SGC_ANI/tmp/*.tmp
		rm ${output_dir}/SGC_ANI/SGC_ANI.txt
	fi

## step-7: tree
	if [ ! -f ${output_dir}/SGC_ANI/UPGMA.nwk ]; then
		log "RUNNING: tree"
		Rscript scripts/tree.r \
		-i ${output_dir}/SGC_ANI/SGC_ANI.tab \
		-o ${output_dir}/SGC_ANI/UPGMA.nwk

		mv Rplots.pdf ${output_dir}/SGC_ANI/

		python scripts/export_phylogenomic_tip_labels.py \
		-i ${output_dir}/SGC_ANI/UPGMA.nwk \
		-o ${output_dir}/SGC_ANI/UPGMA.labels
	fi

## step-8: blast
    mkdir ${output_dir}/blast/files > /dev/null 2>&1
    mapfile -t ids < ${output_dir}/SGC_ANI/UPGMA.labels
    for ((i=0; i<${#ids[@]}-1; i++)); do
        query=${ids[i]}
        subject=${ids[i+1]}

        if [ ! -f ${output_dir}/blast/"${subject}_db.ndb" ]; then
            log "Creating BLAST database for: $subject"

            makeblastdb \
            -in ${output_dir}/cds/${subject}.cds.fasta \
            -dbtype nucl \
            -out ${output_dir}/blast/"${subject}_db"
        fi


        if [ ! -f ${output_dir}/blast/files/"${query}_vs_${subject}.blast" ]; then
            log "Running BLAST: $query vs $subject"

            blastn \
            -query ${output_dir}/cds/${query}.cds.fasta \
            -db ${output_dir}/blast/"${subject}_db" \
            -out ${output_dir}/blast/files/"${query}_vs_${subject}.blast" \
            -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads $cpus
        fi
    done
## step-9: genePlot
    # if [ ! -f ${output_dir}/genePlot.coordinates.png ]; then
        # log "Creating genePlot coordinates"
        # Rscript scripts/gggenomes_plot.coordinates.R \
        # -i ${data_dir}/ \
        # -b ${output_dir}/blast/files/ \
        # -l ${output_dir}/SGC_ANI/UPGMA.labels \
        # --distance 1.0 \
        # --output_format png \
        # --labels_on_top \
        # --char 10 \
        # -o ${output_dir}/genePlot.coordinates.png
    # fi

    # if [ ! -f ${output_dir}/genePlot.genes.png ]; then
        log "Creating genePlot genes"
        Rscript scripts/gggenomes_plot.genes.R \
        -i ${data_dir}/ \
        -b ${output_dir}/blast/files/ \
        -l ${output_dir}/SGC_ANI/UPGMA.labels \
        --distance 0.5 \
        --output_format png \
        --labels_on_top \
        --use_product \
        -o ${output_dir}/genePlot.genes.png
    # fi