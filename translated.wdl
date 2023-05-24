version 1.0

workflow amr_analysis {
    input {
        Array[File] fastqFiles
        Array[Array[File]] pairedReads
        File resfinderDB
    }

    scatter (fastqFile in fastqFiles) {
       call fastqc as fastqcRaw {
           input: 
               fastqFile = fastqFile
       }
    }
    
    call multiqc as multiqcRaw {
       input: 
           input_files = fastqcRaw.outZip
    }
    
    scatter (pairRead in pairedReads) {
        call cutadapt {
            input: read1 = pairRead[0], read2 = pairRead[1]
        }
    }
    

    scatter (fastqFileTrimmed in flatten([cutadapt.outFwd, cutadapt.outRev])) {
       call fastqc as fastqcTrim {
           input: 
               fastqFile = fastqFileTrimmed
       }
    }
    Array[File] fastqcTrimzip = fastqcTrim.outZip
    
    call multiqc as multiqcTrim {
       input: 
           input_files = fastqcTrimzip
    }
    
    # call buildDatabase as resFinder{
    #    input:
    #        database = resfinderDB,
    #        referenceName = "resfinder"
           
    # }

    # scatter (read1 in cutadapt.outFwd){
    #     scatter (read2 in cutadapt.outRev){
    #         call Bowtie2 as resfinderBowtie2{
    #             input:
    #                 R1 = read1,
    #                 R2 = read2,
    #                 database = resFinder.outFile,
    #                 referenceName = "resfinder"
    #         }
    #         call filterBam as resfinderfilterBam{
    #             input:
    #                 samFile = resfinderBowtie2.alignment,
    #                 referenceName = "resfinder"
    #         }

    #         call SortAndIndex as resfinderSortAndIndex{
    #             input:
    #                 bamFile = resfinderfilterBam.bamFile,
    #                 referenceName = "resfinder"
    #         }
    #     }
    # }

    # Array[File] bamFiles = flatten(resfinderSortAndIndex.sortedBam)

    # call CombineResults1 as resfindercombineResults1{
    #         input:
    #             sortedBam = bamFiles[0],
    #             referenceName = "resfinder"
    #     }

    # scatter (sample in bamFiles){
    #     call CombineResults2 as resfindercombineResults2{
    #         input:
    #             sortedReads = sample,
    #             referenceName = "resfinder"
    #     }
    #     call AddSampleNames as resfindersampleNames{
    #         input:
    #             sampleCount = resfindercombineResults2.out,
    #             referenceName = "resfinder"
    #     }
    # }
    # Array [File] renamedSamples = resfindersampleNames.renamedSampleCount

    # call Create_ARG_Genemat{
    #     input:
    #         geneNames = resfindercombineResults1.out,
    #         renamedSampleCounts = renamedSamples,
    #         referenceName = "resfinder"
    # }
 
}

task fastqc {

    input {
        File fastqFile
    }

    Int numCores = 8
    String dockerImage = "quay.io/staphb/fastqc:0.11.9"

    command <<<
        set -euo pipefail

        fastqc -o . ~{fastqFile}
    >>>

    output {
        File outHtml = basename(fastqFile, ".fastq.gz") + "_fastqc.html"
        File outZip = basename(fastqFile, ".fastq.gz") + "_fastqc.zip"
    }

    runtime {
        docker: dockerImage
        disks: "local-disk 500 HDD"
        cpu: numCores
        memory: "16 GB"
    }
}

task multiqc {
    input {
        Array[File]     input_files 

        Boolean         force = false
        Boolean         full_names = false
        String?         title
        String?         comment
        String?         file_name
        String          out_dir = "./multiqc-output"
        String?         template
        String?         tag
        String?         ignore_analysis_files
        String?         ignore_sample_names
        File?           sample_names
        Array[String]?  exclude_modules
        Array[String]?  module_to_use
        Boolean         data_dir = false
        Boolean         no_data_dir = false
        String?         output_data_format
        Boolean         zip_data_dir = false
        Boolean         export = false
        Boolean         flat = false
        Boolean         interactive = true
        Boolean         lint = false
        Boolean         pdf = false
        Boolean         megaQC_upload = false # Upload generated report to MegaQC if MegaQC options are found
        File?           config  # directory
        String?         config_yaml

        String          docker = "quay.io/biocontainers/multiqc:1.8--py_2"
    }

    # get the basename in all wdl use the filename specified (sans ".html" extension, if specified)
    String report_filename = if (defined(file_name)) then basename(select_first([file_name]), ".html") else "multiqc"

    command {
        set -ex -o pipefail

        echo "${sep='\n' input_files}" >  input-filename.txt
        echo "" >>  input-filename.txt

        multiqc \
        --file-list  input-filename.txt \
        --dirs \
        --outdir "${out_dir}" \
        ${true="--force" false="" force} \
        ${true="--fullnames" false="" full_names} \
        ${"--title " + title} \
        ${"--comment " + comment} \
        ${"--filename " + file_name} \
        ${"--template " + template} \
        ${"--tag " + tag} \
        ${"--ignore " + ignore_analysis_files} \
        ${"--ignore-samples" + ignore_sample_names} \
        ${"--sample-names " + sample_names} \
        ${true="--exclude " false="" defined(exclude_modules)}${sep=' --exclude ' select_first([exclude_modules,[]])} \
        ${true="--module " false="" defined(module_to_use)}${sep=' --module ' select_first([module_to_use,[]])} \
        ${true="--data-dir" false="" data_dir} \
        ${true="--no-data-dir" false="" no_data_dir} \
        ${"--data-format " + output_data_format} \
        ${true="--zip-data-dir" false="" zip_data_dir} \
        ${true="--export" false="" export} \
        ${true="--flat" false="" flat} \
        ${true="--interactive" false="" interactive} \
        ${true="--lint" false="" lint} \
        ${true="--pdf" false="" pdf} \
        ${false="--no-megaqc-upload" true="" megaQC_upload} \
        ${"--config " + config} \
        ${"--cl-config " + config_yaml }

        if [ -z "${file_name}" ]; then
            mv "${out_dir}/${report_filename}_report.html" "${out_dir}/${report_filename}.html"
        fi

        tar -c "${out_dir}/${report_filename}_data" | gzip -c > "${report_filename}_data.tar.gz"
    }

    output {
        File multiqc_report            = "${out_dir}/${report_filename}.html"
        File multiqc_data_dir_tarball  = "${report_filename}_data.tar.gz"
    }

    runtime {
        memory: "16 GB"
        cpu: 8
        docker: "${docker}"
        disks: "local-disk 500 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

task cutadapt {
    input {
        File read1
        File read2
    }

    Int numCores = 8
    String dockerImage = "pegi3s/cutadapt"
    String outfw = basename(read1, "001.fastq.gz") + "trimmed.fastq.gz"
    String outrv = basename(read2, "001.fastq.gz") + "trimmed.fastq.gz"

    command {
        set -euo pipefail
        
        cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -O 10 -m 30 -q 20 ~{read1} ~{read2} -o ~{outfw} -p ~{outrv}
    }

    output {
        File outFwd = outfw
        File outRev = outrv
        File outlog = stdout()
    }

    runtime {
        docker: dockerImage
        disks: "local-disk 500 HDD"
        cpu: numCores
        memory: "16 GB"
    }
}

task buildDatabase{
    input{
        File database
        String referenceName
    }
   
    command{
        bowtie2-build ~{database} ${referenceName}
    }
    output {
        Array[File] outFile = glob("~{referenceName}.*.bt2")
    }

    runtime {
        docker: "quay.io/biocontainers/bowtie2:2.4.2--py38h1c8e9b9_1"
    }
}

task Bowtie2 {
    input {
        File R1
        File R2
        Array[File] database
        String referenceName
    }
     String outFile = "${referenceName}_alignment.sam"

    command {
        bowtie2 -x ${referenceName} -1 ~{R1} -2 ~{R2} \
        -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
        --threads 6 | samtools view -Sb -o ~{outFile}
        
    }

    output {
        File alignment = outFile
    }

    runtime {
        docker: "ummidock/innuca"
    }
}

task filterBam {
    input {
        File samFile
        String referenceName
    }

    String outFile = "${referenceName} + output.bam"

    command <<<
      samtools view -h ~{samFile} | awk 'BEGIN {FS="\t"; OFS="\t"} \
    {if (/^@/ && substr($2, 3, 1)==":") {print} \
    else if (($7!="=" || $7=="=") && and($2, 0x40)) {print}}' \
    | samtools view -Shu -o ~{outFile}
    >>>

    output {
        File bamFile = "~{outFile}"
    }

    runtime {
        docker:'ghcr.io/stjudecloud/samtools:1.0.2'
    }
}

task SortAndIndex {
    input {
        File bamFile
        File referenceName
    }

    String outFile = "${referenceName}+sorted.bam"
    String outBai = "${referenceName}+sorted.bam.bai"

    command {
       samtools sort -O bam ~{bamFile} -o ~{outFile}
        samtools index outFile -o~${outBai}
    }

    output {
        File sortedBam = "~{outFile}"
        File bamIndex = "~{outBai}"
    }

    runtime {
        docker: 'ghcr.io/stjudecloud/samtools:1.0.2'
    }
}

task CombineResults1{
    input{
        File sortedBam
        String referenceName
    }
    String outFile = "${referenceName}_out/gene_names.txt"
    command{
        samtools idxstats ~{sortedBam} | grep -v "\*" | cut -f1 -o ~{outFile}
		sed -i '1 i\GENE' ~{outFile}
    }
    output{
        File out = outFile
    }
    runtime {
        docker: 'ghcr.io/stjudecloud/samtools:1.0.2'
    }
}

task CombineResults2{
    input{
        File sortedReads
        String referenceName
    }
    String outFile =  "${referenceName}_out/{sortedReads}_counts.txt"
    command{
        samtools idxstats ${sortedReads} | grep -v "\*" | cut -f3 -o ~{outFile}
    }
    output{
        File out = outFile
    }
    runtime {
        docker: 'ghcr.io/stjudecloud/samtools:1.0.2'
    }
}


task AddSampleNames {
    input {
        File sampleCount
        String referenceName
    }
    String baseName = basename(sampleCount, "_counts")
    String outFile = "renamed_${baseName}_${referenceName}_counts.txt"

    command {
        echo ${baseName} | cat - ${sampleCount} -o ~{outFile}
    }

    output {
        File renamedSampleCount = outFile
    }

}

task Create_ARG_Genemat {
    input {
        File geneNames
        Array[File] renamedSampleCounts
        String referenceName
    }

    File outFile = "${referenceName}_ARG_genemat.txt"

    command <<< 
        paste ${geneNames} ${sep=" " renamedSampleCounts} -o ARG_genemat.txt
    >>>

    output {
        File ARG_genemat = outFile
    }
}