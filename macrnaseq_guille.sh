#!/bin/bash

# Script for RNA-seq processing using nf-core pipeline
# Required software: curl, wget, nextflow, tput, screen

# Parse parameters
while getopts 'n:d:o:x:r:f:g:s:p:' OPTION; do
	case "$OPTION" in
		n)
			# Keep this the same so that of results are found in outdir/project_name/output_per_sample/
			project_name="$OPTARG"
			echo "Project name: $project_name"
			;; # double semicolon to indicate option code has finished
		d)
			# Source directory with all your FASTQ/FQ.gz samples. Add '/*' if you have sample sub-dirs as input_dur/raw_data/*
			input_raw_data_dir="$OPTARG"
			echo "Raw data dir: $input_raw_data_dir"
			;;
		o)
			outdir="$OPTARG"
			echo "Output dir: $outdir"
			;;
		x)
			workdir="$OPTARG"
			echo "Work dir: $workdir"
			;;
		r)
			# If using a non-human/custom reference, add full file paths :
			custom_reference="$OPTARG"
			if [ "$custom_reference" = "y" ] || [ "$custom_reference" = "Y" ]; then
				echo "Using a custom reference genome"
			elif [ "$custom_reference" = "n" ] || [ "$custom_reference" = "N" ]; then
				echo "Using the latest reference genome from Ensembl"
			else
				echo "You must choose 'y' or 'n' after the '-r' flag"
				exit 1
			fi
			;;

		f)
			if [ "$custom_reference" = "y" ] || [ "$custom_reference" = "Y" ]; then
				fasta_file="$OPTARG"
				fasta_filename=$(basename $fasta_file)
				echo "Using custom reference genome: $fasta_filename"
			fi
			;;
		g)
			if [ "$custom_reference" = "y" ] || [ "$custom_reference" = "Y" ]; then
				gtf_file="$OPTARG"
				gtf_filename=$(basename $gtf_file)
				echo "Using custom reference genome annotation file: $gtf_filename"
			fi
			;;
		s)
			if [ $OPTARG = "Hs" ] || [ $OPTARG = "hs" ] || [ $OPTARG = "Human" ] || [ $OPTARG = "human" ] || [ $OPTARG = "Homo_sapiens" ] || [ $OPTARG = "homo_sapiens" ]; then
				species="human"
				echo "Selected species: Homo sapiens"
			elif [ $OPTARG = "Mm" ] || [ $OPTARG = "mm" ] || [ $OPTARG = "Mouse" ] || [ $OPTARG = "mouse" ] || [ $OPTARG = "Mus_musculus" ] || [ $OPTARG = "mus_musculus" ]; then
				species="mouse"
				echo "Selected species: Mus musculus"
			else
				echo "Supported species are human (Hs, hs, Human, human, Homo_sapiens or homo_sapiens) \\
				and mouse (Mm, mm, Mouse, mouse, Mus_musculus or mus_musculus)."
				echo "script usage: $basename \$0) [-n project_name] [-d raw_data_dir] [-o outdir] [-x workdir] \\
	                        [-r use_custom_ref_genome] [-f fasta_file] [-g gtf_file] -s [species]" >&2 # ">&2 redirects normal output to stderr"
				exit 1
			fi
			;;
		p)
			# File with parameters for nf-core/rnaseq, optional
			params_file="$OPTARG"
			echo "Using parameters file: $params_file"
			;;
		?)
			echo "script usage: $basename \$0) [-n project_name] [-d raw_data_dir] [-o outdir] [-x workdir] \\
			[-r use_custom_ref_genome] [-f fasta_file] [-g gtf_file] -s [species]" >&2 # ">&2 redirects normal output to stderr"
			exit 1
			;;
	esac
done
shift "$(($OPTIND -1))"


# Exit on error, undefined variable, or a pipe failure
set -euo pipefail

################################################################################

# Error checking
check_inputs () {

	# test -z "$var"  returns true if length of variable == 0
	if test -z "$project_name"; then 
		echo "Error: 'project_name' cannot be empty. Please provide a valid project name."
				exit
			fi
		
	if test -z "$input_raw_data_dir"; then 
		echo "Error: 'input_raw_data_dir' cannot be empty. Please provide the source directory for FASTQ/FQ.gz samples."
		exit
	fi

	if test -z "$outdir"; then	
	echo "Error: 'outdir' cannot be empty. Please provide an output directory."
		exit
	fi

	if test -z "$workdir"; then 
		echo "Error: 'workdir' cannot be empty. Please provide a working directory."
	exit
fi

}

# Function to check if a file exists
file_exists () {
	[[ -f "$1" ]]
}

# Download reference files
download_reference_files() {

	echo "Downloading reference files... (~1000 MB)"
echo $species
	
	if [ $species = "human" ]; then
		echo "Downloading latest H. sapiens reference genome FASTA and GTF files..."
	wget -q --show-progress -c -L "ftp://ftp.ensembl.org/pub/release-${latest_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz" -O "$fasta_file"
	wget -q --show-progress -c -L "ftp://ftp.ensembl.org/pub/release-${latest_release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${latest_release}.gtf.gz" -O "$gtf_file"
elif [ $species = "mouse" ]; then
		echo "Downloading latest M. musculus reference genome FASTA and GTF files..."
		wget -q --show-progress -c -L "ftp://ftp.ensembl.org/pub/release-${latest_release}/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz" -O "$fasta_file"
		wget -q --show-progress -c -L "ftp://ftp.ensembl.org/pub/release-${latest_release}/gtf/mus_musculus/Mus_musculus.GRCm39.${latest_release}.gtf.gz" -O "$gtf_file"
	fi
	
	echo
	echo "$(tput setaf 5)All files succesffuly dowloaded and stored in:$(tput sgr0) $download_dir/reference"
}

# Function to check reference file versioning and download if new needed
fetch_reference_files () {
    echo "Fetching reference files..."

    # Check if a custom reference is provided
    if [[ "$custom_reference" == "y" ]]; then
        echo "$(tput setaf 3)Custom reference files selected.$(tput sgr0)"
        if [[ -n "$fasta_file" && -n "$gtf_file" ]]; then
            if file_exists "$fasta_file" && file_exists "$gtf_file"; then
                echo "$(tput setaf 69)Custom GTF/FASTA reference files provided and exist.$(tput sgr0)"
            else
                echo "$(tput setaf 1)Error: Custom reference files not found. Please check the file paths.$(tput sgr0)"
                exit 1
            fi
        else
            echo "$(tput setaf 1)Error: Custom reference file paths are not provided.$(tput sgr0)"
            exit 1
        fi
    else
        # Proceed with fetching reference files
        # Fetch the latest Ensembl release number
        latest_release=$(curl -s 'http://rest.ensembl.org/info/software?content-type=application/json' | grep -o '"release":[0-9]*' | cut -d: -f2)

        # File paths for reference files
        download_dir=$PWD
        mkdir -p "$download_dir/reference"

        if [ $species = "human" ]; then
        fasta_file="$download_dir/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
        gtf_file="$download_dir/reference/Homo_sapiens.GRCh38.${latest_release}.gtf.gz"
elif [ $species = "mouse" ]; then
		fasta_file="$download_dir/reference/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"
		gtf_file="$download_dir/reference/Mus_musculus.GRCm39.${latest_release}.gtf.gz"
	fi

        # Check if reference files exist or download them
        if file_exists "$fasta_file" && file_exists "$gtf_file"; then
            echo "$(tput setaf 69)GTF/FASTA reference files are already available.$(tput sgr0)"
        else
            echo ">> New versions of the reference genome GTF and FASTA are available."
			echo ">> Reference files are required to run nf-core/rnaseq."
			read -p "	>> Do you want to download them? (y/n): " choice
            case "$choice" in 
                y|Y ) download_reference_files ;;
                n|N ) echo "Files not downloaded. $(tput setaf 3)The older version reference genome will be used.$(tput sgr0)";;
                * ) echo "Invalid choice. Files not downloaded." ;;
            esac
        fi
    fi
}

# Check if user is in a screen session
linux_screen_check () {
	echo "Setting up background envirnoment..."
	while true; do
		read -p "Are you in a Linux screen session? (y/n): " response
		case "$response" in
			[yY]* ) break ;;
			[nN]* ) 
				echo
				echo "You $(tput setaf 1)MUST$(tput sgr0) set up a Linux screen, otherwise the run will halt."
				echo
				echo "$(tput smul)Quick guide:$(tput sgr0)"
				echo "Use $(tput setaf 4)screen -S screen_name$(tput sgr0) to create a new screen."
				echo "Use $(tput setaf 4)screen -r screen_name$(tput sgr0) to re-attach to an existing screen."
				echo "	When attached to a screen:"
				echo "		use Ctrl+A+D to detach it"
				echo "		use Ctrl+A+Esc to scroll and then press q to stop"
				echo "Use $(tput setaf 4)screen -X -S screen_name quit$(tput sgr0) to delete an existing screen."
				echo
				echo "Create/re-attach to a linux screen and re-run to carry on"
				echo 
				echo "                         B Y E"
				echo
				exit 1 ;;
			* ) echo "Invalid response. Please enter 'y' for yes or 'n' for no.";;
		esac
	done
	echo "$(tput setaf 69)Background environment is set up$(tput sgr0)"
}

# Prepare samplesheet
prepare_samplesheet () {
        echo "Preparing samplesheet..."
        echo "Raw data source: $(tput setaf 3)$input_raw_data_dir$(tput sgr0)"

        # Create array with sample names
        sample_names=( $(find $input_raw_data_dir/*q.gz -type f -exec basename {} ';' | cut -d "_" -f 1 | uniq) )

        # List fastq files containing R1 and R2
        read1=( $(ls $input_raw_data_dir/*R1*q.gz) )
        read2=( $(ls  $input_raw_data_dir/*R2*q.gz) )

	nfiles_read1=${#read1[@]} # ${#array[@]} gets the array length
	nfiles_read2=${#read2[@]}
	
	echo "$nfiles_read1 files found for R1"
	echo "$nfiles_read2 files found for R2"
	
	if [ $nfiles_read1 != $nfiles_read2 ]; then
		echo "ERROR: The number of files for R1 and R2 is different"
		exit 1
	fi

        n_samples=${#sample_names[@]}
        echo "Uniquely named samples found on raw data source >> $(tput setaf 2)$n_samples$(tput sgr0) <<"

        # Build sample sheet
        echo "sample,fastq_1,fastq_2,strandedness" > samplesheet.csv # Header line
        for ((i=0;i<$n_samples;i++)); do echo -e ${sample_names[i]}","${read1[i]}","${read2[i]}",auto" >> samplesheet.csv; done # rest of lines

        samplesheet_nsamples=$(tail -n+2 samplesheet.csv | cut -d , -f 1 | wc -l )
        echo "Uniquely named samples found on samplesheet.csv >> $(tput setaf 2)$samplesheet_nsamples$(tput sgr0) <<"

        if [[ $samplesheet_nsamples != $n_samples ]]; then
                echo "Not all samples were correctly added to the sample sheet!"
                else
                echo "$(tput setaf 69)Samplesheet made successfully$(tput sgr0)"
        fi
}


# Nf-core/rnaseq parameters. Adjust as needed
nf_core_rnaseq_command () {

	# If no configuration file was provided, run with default values
	if test -z "$params_file"; then
	
		echo "$(tput setaf 69)Using script default paramters...$(tput sgr0)"
	
		nextflow run nf-core/rnaseq \
			-profile docker \
			--genome null \
			--igenomes_ignore \
			--gtf "$gtf_file" \
			--fasta "$fasta_file" \
			--input input.csv \
			--outdir "$sample_outdir" \
			-w "$sample_workdir" \
			-with-tower
		nextflow clean -f -q
	elif file_exists "$params_file"; then
	
		echo "$(tput setaf 69)Using custom configuration file...$(tput sgr0)"
	
		nextflow run nf-core/rnaseq \
			-profile docker \
			-params-file "$params_file" \
			--genome null \
			--igenomes_ignore \
			--gtf "$gtf_file" \
			--fasta "$fasta_file" \
			--input input.csv \
			--outdir "$sample_outdir" \
			-w "$sample_workdir" \
			-with-tower
			
		nextflow clean -f -q
	fi
}

# Run RNA-Seq on samples
run_nf_core_pipeline () {

	echo "Running nf-core/rnaseq..."

	input="samplesheet.csv"
	processed_samples_log="processed_samples.log"
	touch "$processed_samples_log"

	head -1 "$input" > head.tmp
	tail -n+2 "$input" | cut -d , -f 1 | sort | uniq > samples.tmp
	samples_to_process_count=$(wc -l < samples.tmp)

	
	i=1
	for sample in $(cat samples.tmp); do

		# Set up working and output directories
		sample_workdir="$workdir/$project_name/$sample"
		sample_outdir="$outdir/$project_name/$sample"
		mkdir -p "$sample_outdir"
		mkdir -p "$sample_workdir"


		if grep -q "Workflow execution completed successfully!" $sample_outdir/pipeline_info/execution_report_*.html 2>/dev/null
		then
			echo "$(tput setaf 2)The sample $sample excution completed successfuly $(tput sgr0)"
			((samples_to_process_count--))
			continue
		fi

		echo "$(tput setaf 3)Processing $sample... ($i/$samples_to_process_count)$(tput sgr0)"

		# Prepare input.csv for each sample
		cat head.tmp > input.csv
		awk -F',' -v sample="$sample" '$1 == sample' samplesheet.csv >> input.csv

		# Run the nf-core/rnaseq pipeline
		nf_core_rnaseq_command

		rm input.csv
		((i++))
	done
	
	# Clean up
	rm head.tmp samples.tmp 
	echo "$(tput setaf 69)RNA-Seq processing of all $project_name samples completed successfully $(tput sgr0)"
}

out_message () {

# Cool message using echo
	echo "$(tput setaf 2)THANK YOU$(tput sgr0) for using the $(tput setaf 69)MAC-RNASEQ$(tput sgr0) pipeline!"
	echo "$(tput sitm)Hope it was easy to work with :)$(tput sgr0)"
	echo "$(tput bold)Please$(tput sgr0) delete $(tput setaf 3)$workdir$(tput sgr0) after all further analysis is completed!"
}

################################################################################

# Begin script
# Comment out unncessary steps as needed

echo "---------------------------------------------------- SET UP ---"
echo "$(tput setaf 90)MAC-RNASEQ$(tput sgr0)"
		check_inputs
echo "---------------------------------------------------- STEP 1 ---"
		fetch_reference_files
echo "---------------------------------------------------- STEP 2 ---"
		linux_screen_check
echo "---------------------------------------------------- STEP 3 ---"
		prepare_samplesheet
echo "---------------------------------------------------- STEP 4 ---"
		# Comment the pipeline line out to try out other features first
		#run_nf_core_pipeline
echo "---------------------------------------------------------------"
		out_message
