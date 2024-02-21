#!/usr/bin/env nextflow
@Grab('com.xlson.groovycsv:groovycsv:1.1')
import static com.xlson.groovycsv.CsvParser.parseCsv
import java.io.File;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                     Checkup imported Variables                                                          //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//Variables will be needed for configuration of the programm, maybe a switch to configuration file coming soon  
//Folders of interest, annotation and alignment files
println "General Threads: "
//params.threads = 30
println params.threads
println " "

//params.general_folder = "/media/stefan/Samsung_T5/test_data/"
println "General folder: "
println params.general_folder
println " "

params.data_folder = []
params.sample_names = []
params.suffix

//println "Mapper"
//println params.mapper

//Define a path that contains the file you want to align your data to. Usually this file is in .fasta format
//params.genome_fasta = "/media/stefan/Samsung_T5/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
println "Genome Fasta chosen: "
println params.genome_fasta 
println " "

//Define a path that contains the transcriptome alignment file
println "Transcriptome Fasta chosen: "
println params.transcriptome_fasta 
println " "


//Define a path that contains the file you want to annotate your data to. Usually this file is in .gtf format
//params.genome_gtf = "/media/stefan/Samsung_T5/Homo_sapiens.GRCh38.104.gtf" 
println "Gtf File chosen: "
println params.genome_gtf
println " "

//Define a directory in which all the scripts belonging to the application are stored. Especially fc.py and merge_fc.py are important for this script.
//params.script_dir = "/home/stefan/arbeit_gerber/nextflow_dRNA_seq/"
println "Script directory: "
println params.script_dir
println " "


//Define a directory in which you want to memorize all your run specific data and the progress status of your data. This directory will be automatically created if it does not exist so far
//params.run_dir = "/media/stefan/Samsung_T5/run_dir2/"
println "Run directory: "
println params.run_dir
println " "

//Conditions that will be compared in DESeq
println "Conditions compared: "
println params.conditions
println " "

//Are samples barcoded
println "Barcoded: "
println params.barcoded
println " "

//Metadata file
println "Metadata: "
println params.metadata
println " "

println "BED File chosen: "
println params.bed_file
println " "

params.batchsize = 30
println "Maximal batch size: "
println params.batchsize
println " "
//Channels are important for Activating and reactivating Processes at a given timepoint

//params.DRS = 1
println "Type of Sample:"
println params.DRS
println " "

checkup = Channel.create()
//checkup2 = Channel.create()
//checkup3 = Channel.create()
done = Channel.create()
done2 = Channel.create()


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                       Reference Variables                                                               //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Definition of a global counter class to ensure copy by reference operations on counter

class IntObj {
    public int value;
}

counter = new IntObj()
counter.value 

counter2 = new IntObj()
counter2.value = 0

counter_iterations = new IntObj()
counter_iterations.value = 0


alignment_running_counter = new IntObj()
alignment_running_counter.value = 0

alignment_running_counter_transcript = new IntObj()
alignment_running_counter_transcript.value = 0

fc_running = new IntObj()
fc_running.value = 0

suffix_determined = new IntObj()
suffix_determined.value = 0

fc_merging = new IntObj()
fc_merging.value = 0

deseq_running = new IntObj()
deseq_running.value = 0

move_transcripts_minimap2_running = new IntObj()
move_transcripts_minimap2_running.value = 0

salmon_annotation_running = new IntObj()
salmon_annotation_running.value = 0

continue_bool = new IntObj()
continue_bool.value = 0

iteration = new IntObj()
iteration.value = 0

//Process finds out which sample folders exist
process find_data_folders(){
    input: 

    output:
    val done into data_folder_found
    exec:
    if (params.barcoded == 1){
        def text = new File(params.metadata).getText()
        def data = parseCsv(text, separator: '\t')
        data.each{line -> params.sample_names.add(line["Samples"])}
        //println params.sample_names
    }
    else {
        def text = new File(params.metadata).getText()
        def data = parseCsv(text, separator: '\t')
        data.each{line -> params.sample_names.add(line["Samples"])}
        //println params.sample_names
    }
    done = 1
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                            Make crucial Directories                                                                     //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Process creates necessary directories in each folder of interest respectively, performs minimap indexing, and performs a structure conversion of the gtf to a csv file
process make_directories{
    input:
    val done2 from data_folder_found
    output:
	val finished into start_process_channel
    val finished into side_processing1
    val finished into side_processing2
	
	script:
	string = ""
	for(i in 0..params.sample_names.size()-1){
	    string = string + params.run_dir + params.sample_names.get(i) + "/ "
	}
	println "Make directories in ${string}"
	finished = 1
	"""
    if [ ! -d ${params.run_dir} ]
	then
	    mkdir ${params.run_dir}
	fi

    if [ ! -d ${params.run_dir}r_objects ]
    then   
        mkdir ${params.run_dir}r_objects
    fi

    if [ ! -d ${params.run_dir}error_logs ]
    then
        mkdir ${params.run_dir}error_logs
    fi

    if [ ! -d ${params.run_dir}bam_genome_merged ]
    then
        mkdir ${params.run_dir}bam_genome_merged
    fi
    if [ ! -d ${params.run_dir}bam_transcriptome_merged ]
    then
        mkdir ${params.run_dir}bam_transcriptome_merged
    fi
    if [ ! -d ${params.run_dir}ReadLengthFolder ]
    then
        mkdir ${params.run_dir}ReadLengthFolder
    fi
	for i in ${string}
	do
      if [ ! -d \${i} ]
      then
         mkdir \${i}
      fi

      if [ ! -d \${i}bam_files ]
      then
        mkdir \${i}bam_files
      fi

      if [ ! -d \${i}bam_files_transcripts ]
      then
        mkdir \${i}bam_files_transcripts 
      fi

      if [ ! -d \${i}bam_files_transcripts/full ]
      then
            mkdir \${i}bam_files_transcripts/full
      fi

	  if [ ! -d \${i}merged_fastq ]
	  then
	     mkdir \${i}merged_fastq
	  fi 

	  if [ ! -d \${i}merged_fc ]
	  then
	     mkdir \${i}merged_fc
	  fi 

      if [ ! -d \${i}merged_fc_splice ]
	  then
	     mkdir \${i}merged_fc_splice
	  fi 

      if [ ! -d \${i}single_fc ]
      then 
         mkdir \${i}single_fc
      fi

      if [ ! -d \${i}single_fc_splice ]
      then 
         mkdir \${i}single_fc_splice
      fi
      
      if [ ! -d \${i}fastqc ]
      then
         mkdir \${i}fastqc
      fi
      
      if [ ! -d \${i}salmon ]
      then
         mkdir \${i}salmon
      fi

	done
    echo 1 > ${params.run_dir}process_running.txt

    if [ ${params.DRS} -eq 1 ]
    then
        minimap2 --MD -ax splice -uf -k14 -d ${params.run_dir}MT-human_ont.mmi $params.genome_fasta || echo "Indiexing genome failed" >> ${params.run_dir}error_logs/index_genome.log
        minimap2 --MD -ax map-ont -uf -k14 -d ${params.run_dir}MT-human_transcript_ont.mmi $params.transcriptome_fasta || echo "Indexing transcriptome failed" >> ${params.run_dir}error_logs/index_transcriptome.log
    else
        minimap2 --MD -ax splice -d ${params.run_dir}MT-human_ont.mmi $params.genome_fasta || echo "Indiexing genome failed" >> ${params.run_dir}error_logs/index_genome.log
        minimap2 --MD -ax map-ont -d ${params.run_dir}MT-human_transcript_ont.mmi $params.transcriptome_fasta || echo "Indexing transcriptome failed" >> ${params.run_dir}error_logs/index_transcriptome.log
    fi
    """
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                Feature percentiles and gtf conversion                                                   //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process side_processes{
    input:
    val finished from side_processing1

    script:
    """
    python ${params.script_dir}convert_gtf_to_df.py -i ${params.genome_gtf} -o ${params.run_dir}converted_gtf.csv || echo "Conversion failed" >> ${params.run_dir}error_logs/gtf_conversion.log
    """
}

process side_processes2{
    input:
    val finished from side_processing2

    script:
    """
    python ${params.script_dir}createFeaturePercentiles.py --bed ${params.bed_file} --output_dir ${params.run_dir} || echo "Feature percentile preprocess not working" >> ${params.run_dir}error_logs/feature_percentile_preprocess_failed.log
    """
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                Progress lookup in run folder                                                            //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process see_progress{
	input:
	val finished from start_process_channel
	
	output:
	val done into activate_channels

	script:
	File file = new File("${params.run_dir}data_seen.txt")
	if(file.exists()){
    		    println "Data seen exists"
    		    def line = 0;
    		    params.data_seen_list = []
                file.withReader{reader -> while((line = reader.readLine()) != null){
     		    data = line.split(" ").collect{it as String}
                if (data.size() > 0){
                    for (i in 0..data.size()-1){
                        params.data_seen_list.add(data.get(i))
                        }
                    }
                }
                }
     		    //println "Data seen list: ${params.data_seen_list} "
                if (params.data_seen_list.size() == 0){
                    file.write(" ")
                    params.data_seen_list = []
                }
	    
    }
	else{
   	    println "Data seen does not exist"
   	    params.data_seen_list = []
   	    file.write(" ")
	}
    File file2 = new File("${params.run_dir}current_iteration.txt")
    if(file2.exists()){
    		    def line = 0;
                file2.withReader{reader -> while((line = reader.readLine()) != null){
     		    data = line.split(" ").collect{it as int}
                if (data.size() > 0){
                    iteration.value = data.get(0)
                }
                }
        }
    }
    else
    {
        iteration.value = 0
    }
    println "Current iteration $iteration.value"
    """
    iteration=${iteration.value}
    if [ \$iteration -gt 0 ]
    then
    iteration=\$((\${iteration}+1))
    awk -F "," -v iter=\$iteration '{ if (\$2 < iter) print \$1,\$2,\$3 }' ${params.run_dir}processing_time_table.csv | sed 's/ /,/g;w ${params.run_dir}processing_time_table2.csv'
    echo Tool,Iteration,Time > ${params.run_dir}processing_time_table3.csv & cat ${params.run_dir}processing_time_table2.csv >> ${params.run_dir}processing_time_table3.csv
    cp ${params.run_dir}processing_time_table3.csv ${params.run_dir}processing_time_table.csv
    rm ${params.run_dir}processing_time_table2.csv
    rm ${params.run_dir}processing_time_table3.csv
    fi
    """


}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                 Starting the analysis loop                                                              //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process activate_channels{
     input:
     val done from activate_channels
     
     exec:
     //Activation of first processes by adding a value into a channel
     checkup << Channel.from(1)
     //checkup2 << Channel.from(1)

}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                        Detection of data in paths of interest                                                           //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





//Process checks for all files in each folder of interest respectively 

process fromPath_data{
    memory "4 GB"
    input:
    val check from checkup

    output:
    val files into alignment_prep


    exec:
    
    //While loop is used to prevent from major RAM use in case of to many parallel alignment processes
	println "Alignment running"
    while(alignment_running_counter.value == 1){
	    sleep(3000)
	}
    println "Extracting files from Path"
    files = []
    if (params.barcoded == 1){    
        try{
            files << Channel.fromPath("${params.general_folder}*/*/fastq_pass/*/*", type: "file").collect()
        }
        catch(Exception e){}
    }
    else{
        try{ 
            files << Channel.fromPath("${params.general_folder}*/*/fastq_pass/*", type: "file").collect()
            //println files
        }
        catch(Exception e){}  
    }
        //println files
    done << Channel.from("yes")
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                       Preprocessing                                                                     //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Process compares value of data_seen_list with the list of all files that are present in the folders of interest at each iteration and extracts new files. An array of new files and a concatenated string 
//with format "file1 data2 data3" is transported to the next processes

process data_alignment_prep{
        memory "4 GB"
        input:
        val files from alignment_prep
        
        output:
        val string into minimap_channel
        val string_array into minimap_channel2
        val string into minimap_transcript_channel
        val string_array into minimap_transcript_channel2
        val string into merge_channel
        val string into size_checkup
        val string_array into size_checkup2
        val data_string into merge_fastq_channel
        val data_seen_string into merge_fastq_channel2
        val string into merge_fastq_channel3
        

        script:

        string = ""
        string_array = []
        println "Preprocessing \n"
        println "Data seen before update: ${params.data_seen_list.size()} \n"
        println ""
        data_to_process_list = []
        if (files.val.get(0).getClass() == nextflow.util.ArrayBag){
        for(j in 0..files.val.size()-1){
            for(k in 0..files.val.get(j).size()-1){
                if(files.val.get(j).get(k).toString() in params.data_seen_list){
                }
                else if(files.val.get(j).get(k) in params.data_seen_list){
                }
                else{
                    data_to_process_list.add(files.val.get(j).get(k))
                }
            }
        }
        
        temp_list = []
        if (data_to_process_list.size > 1){
        for (q in 0..data_to_process_list.size()-1){
            if (params.barcoded == 1){    
               if (data_to_process_list.get(q).getParent().name in params.sample_names){
                   temp_list.add(data_to_process_list.get(q))
               }
            }
            else{
               if (data_to_process_list.get(q).getParent().getParent().getParent().name in params.sample_names){
                   temp_list.add(data_to_process_list.get(q))
               }
            }
        }
        data_to_process_list = temp_list
        }
        }
        if (suffix_determined.value == 0){
            if(data_to_process_list.size > 0){
                first_file = data_to_process_list.get(0).toString()
                first_file_parts = first_file.split('\\.').collect()
                last_part = first_file_parts.get(first_file_parts.size()-1)
                println "Last part: "
                println last_part
                println ""
                if (last_part == "gz"){
                    params.suffix = "fastq.gz"
                }
                else {
                    params.suffix = "fastq"
                }
                suffix_determined.value = 1
                println "Params suffix: "
                println params.suffix
                println ""
            }
        }
        println "Data to process: ${data_to_process_list.size()}"
        x_seen = []
        for(k in 0..Math.min(params.batchsize,data_to_process_list.size()-1)){
            if (data_to_process_list.size() > 1){
                x = Math.abs(new Random().nextInt() % (data_to_process_list.size()))
                //While loop is integrated to avoid double insertion of one file 
                while(x in x_seen){
                    x = Math.abs(new Random().nextInt() % (data_to_process_list.size()))
                }
                x_seen.add(x)
                params.data_seen_list.add(data_to_process_list.get(x))
                string = string + data_to_process_list.get(x) + " "
                string_array.add(data_to_process_list.get(x))
            }
            else if (data_to_process_list.size() == 1){
                params.data_seen_list.add(data_to_process_list.get(0))
                string = string + data_to_process_list.get(0) + " "
                string_array.add(data_to_process_list.get(0))   
            }
            else {
                string = ""
                string_array = []
            }
        }
        if (string != ""){
            alignment_running_counter.value = 1
            iteration.value = iteration.value + 1
        }
        println ""
        println "New data: ${string_array.size()} \n"
        println "Data seen: ${params.data_seen_list.size()}"
        data_seen_string = ""
        data_string = ""  
        File file_tmp = new File("${params.run_dir}data_seen_tmp_processed.txt")
        file_tmp.write("")
        for(i in 0..params.data_seen_list.size()-1){
            if (params.data_seen_list.size() > 0){
                if (params.data_seen_list.get(i).toString() != ""){
                    data_seen_string = data_seen_string + params.data_seen_list.get(i).toString() + " "
                    file_tmp.append(params.data_seen_list.get(i).toString() + "\n")
                }
            }
        }
        for(i in 0..params.sample_names.size()-1){
        data_string = data_string + params.sample_names.get(i) + " "
        }
        """
        corruption_gene_bam=0
        corruption_transcript_bam=0
        bam_files_seen=" "
        bam_files_seen_transcripts=" "
        for ((i=1; i<=\$(cat ${params.run_dir}data_seen_tmp_processed.txt | wc -l); i++))
        do
            I=\$(sed -n \$((\$i))p ${params.run_dir}data_seen_tmp_processed.txt)
            if [ ${params.barcoded} -eq 1 ]
            then
                basis=\$(basename \$(dirname \${I}))
                filename=\$(basename \${I})
                converted_filename=\${filename/${params.suffix}/bam}
                outname=${params.run_dir}\${basis}/bam_files/\${converted_filename}
                outname_transcripts=${params.run_dir}\${basis}/bam_files_transcripts/\${converted_filename}
            else
                basis=\$(basename \$(dirname \$(dirname \$(dirname \${I}))))
                filename=\$(basename \${I})
                converted_filename=\${filename/${params.suffix}/bam}
                outname=${params.run_dir}\${basis}/bam_files/\${converted_filename}
                outname_transcripts=${params.run_dir}\${basis}/bam_files_transcripts/\${converted_filename}
            fi
            bam_files_seen=\${bam_files_seen}\${outname}" "
            bam_files_seen_transcripts=\${bam_files_seen_transcripts}\${outname_transcripts}" "
        done
        #echo \$bam_files_seen > ${params.run_dir}error_logs/bamfilesseen_bam_files.txt 
        #echo \$bam_files_seen_transcripts > ${params.run_dir}error_logs/bamfilesseen_transcripts_bam_files.txt 
        for i in ${data_string}
        do
            files_in_i=\$(ls ${params.run_dir}\${i}/bam_files/)
            files_in_i_transcripts=\$(ls ${params.run_dir}\${i}/bam_files_transcripts/)
            #echo \$files_in_i >> ${params.run_dir}/error_logs/filesini_bam_files.txt \${i}
            #echo \$files_in_i >> ${params.run_dir}/error_logs/filesini_transcript_bam_files.txt \${i}
            for j in \$files_in_i
                do
                    if [[ ! "\$bam_files_seen" =~ "${params.run_dir}\${i}/bam_files/\${j}" ]]
                    then
                        if [[ ! -d ${params.run_dir}\${i}/bam_files/\${j} ]]
                        then
                            rm ${params.run_dir}\${i}/bam_files/\${j} && echo ${params.run_dir}\${i}/bam_files/\${j} >> ${params.run_dir}error_logs/removed_bam_files.txt
                            corruption_gene_bam=1
                            echo "1" > ${params.run_dir}corruption_cleanup.txt
                        fi
                    fi
                done
            for j in \$files_in_i_transcripts
                do
                    if [[ ! "\$bam_files_seen_transcripts" =~ "${params.run_dir}\${i}/bam_files_transcripts/\${j}" ]]
                    then
                        if [[ ! -d ${params.run_dir}\${i}/bam_files_transcripts/\${j} ]]
                        then
                            rm ${params.run_dir}\${i}/bam_files_transcripts/\${j} && echo ${params.run_dir}\${i}/bam_files_transcripts/\${j} >> ${params.run_dir}error_logs/removed_bam_files_transcripts.txt
                            corruption_transcript_bam=1
                        fi
                    fi
                done
        done
        echo \$corruption_gene_bam >> ${params.run_dir}error_logs/corruption_gene.log
        function samtoolsParallelAll  {
            if [ ! \$(ls ${params.run_dir}\$1/bam_files/*.bam | wc -l) -eq 0 ]
            then
                par_basis=\$1
                run_dir=\$2
                sample=\${run_dir}\${par_basis}/bam_files/
                echo "Running samtoolsParallelAll" >> ${params.run_dir}error_logs/samtoolsParallelall_running.log
                sample_transcripts=\${run_dir}\${par_basis}/bam_files_transcripts/
                #samtools merge \${run_dir}\${par_basis}/all_gene.bam \${sample}*.bam -f -h \${sample}*.bam --threads ${params.threads} -c -p || echo "\${sample}" >> \${run_dir}error_logs/bam_file_genome_merge_after_corruption_error.log
                samtools merge \${run_dir}\${par_basis}/all_gene.bam \${sample}*.bam -f --threads ${params.threads} -c -p || echo "\${sample}" >> \${run_dir}error_logs/bam_file_genome_merge_after_corruption_error.log
                cp \${run_dir}\${par_basis}/all_gene.bam \${run_dir}bam_genome_merged/\${par_basis}.bam
                samtools index \${run_dir}bam_genome_merged/\${par_basis}.bam
                #samtools merge \${run_dir}\${par_basis}/salmon/all.bam \${sample_transcripts}*.bam -f -h \${sample_transcripts}*.bam --threads ${params.threads} -c -p || echo "\${sample_transcripts}" >> \${run_dir}error_logs/bam_file_transcriptome_merge_after_corruption_error.log
                samtools merge \${run_dir}\${par_basis}/salmon/all.bam \${sample_transcripts}*.bam -f --threads ${params.threads} -c -p || echo "\${sample_transcripts}" >> \${run_dir}error_logs/bam_file_transcriptome_merge_after_corruption_error.log
                cp \${run_dir}\${par_basis}/salmon/all.bam \${run_dir}bam_transcriptome_merged/\${par_basis}.bam
                samtools index \${run_dir}bam_transcriptome_merged/\${par_basis}.bam
            fi
        }
        export -f samtoolsParallelAll
        function fastqmergeParallel {
        par_basis=\$1
        run_dir=\$2
        data_to_process=""
        for ((i=1; i<=\$(cat ${params.run_dir}data_seen_tmp_processed.txt | wc -l); i++))
        do
            I=\$(sed -n \$((\$i))p ${params.run_dir}data_seen_tmp_processed.txt)
            if [ ${params.barcoded} -eq 1 ]
            then
            basis=\$(basename \$(dirname \$I))
            else
            basis=\$(basename \$(dirname \$(dirname \$(dirname \$I))))
            fi
            if [ \${basis} == \${par_basis} ]
            then
                data_to_process=\$data_to_process\$I" "
            fi
        done
        if [ ${params.barcoded} -eq 0 ]
        then
            if [ ! \$(ls ${params.general_folder}\${par_basis}/*/fastq_pass/*.fastq.gz | wc -l) -eq 0 ]
                then
                zcat \${data_to_process} | gzip > \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq.gz
                fi
            if [ ! \$(ls ${params.general_folder}\${par_basis}/*/fastq_pass/*.fastq | wc -l) -eq 0 ]
                then
                cat \${data_to_process} > \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq
                fi
        else
            if [ ! \$(ls ${params.general_folder}*/*/fastq_pass/\${par_basis}/*.fastq.gz | wc -l) -eq 0 ]
            then
            zcat \${data_to_process} | gzip > \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq.gz
            fi
            if [ ! \$(ls ${params.general_folder}*/*/fastq_pass/\${par_basis}/*.fastq | wc -l) -eq 0 ]
            then
            cat \${data_to_process} > \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq
            fi
        fi      
        }
        export -f fastqmergeParallel 
       

        if [ -f ${params.run_dir}corruption_cleanup.txt ]
        then
            cleanup=\$(cat ${params.run_dir}corruption_cleanup.txt)
            if [ \$cleanup -eq 1 ]
            then 
            parallel -v -u --env samtoolsParallelAll --no-notice -j ${params.threads} samtoolsParallelAll ::: ${data_string} ::: ${params.run_dir}
            parallel -v -u --env fastqmergeParallel --no-notice -j ${params.threads} fastqmergeParallel ::: ${data_string} ::: ${params.run_dir}
            fi
        fi
        if [ ! -f ${params.run_dir}processing_time_table.csv ]
            then
                echo "Tool,Iteration,Time" >> ${params.run_dir}processing_time_table.csv
            fi
        iteration=${iteration.value}
        echo \$((\${iteration})) > ${params.run_dir}/current_iteration.txt
        """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                          Merge fastq                                                                    //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


process merge_fastq{
    input:
    val data_string from merge_fastq_channel
    val data_seen_string from merge_fastq_channel2
    val string from merge_fastq_channel3
    output:
    val data_string into fastq_merge_done_channel

    script:
    if (string != "")
    """
    function fastqmergeParallel {
        par_basis=\$1
        run_dir=\$2
        data_to_process=""
        for i in ${string}
        #for ((i=1; i<=\$(cat ${params.run_dir}data_seen_tmp_processed.txt | wc -l); i++))
        do
            I=\$i
            #I=\$(sed -n \$((\$i))p ${params.run_dir}data_seen_tmp_processed.txt)
            if [ ${params.barcoded} -eq 1 ]
            then
            basis=\$(basename \$(dirname \$I))
            else
            basis=\$(basename \$(dirname \$(dirname \$(dirname \$I))))
            fi
            if [ \${basis} == \${par_basis} ]
            then
                data_to_process=\$data_to_process\$I" "
            fi
        done
        if [ ${params.barcoded} -eq 0 ]
        then
            if [ ! \$(ls ${params.general_folder}\${par_basis}/*/fastq_pass/*.fastq.gz | wc -l) -eq 0 ]
                then
                zcat \${data_to_process} | gzip -c >> \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq.gz
                fi
            if [ ! \$(ls ${params.general_folder}\${par_basis}/*/fastq_pass/*.fastq | wc -l) -eq 0 ]
                then
                cat \${data_to_process} >> \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq
                fi
        else
            if [ ! \$(ls ${params.general_folder}*/*/fastq_pass/\${par_basis}/*.fastq.gz | wc -l) -eq 0 ]
            then
            zcat \${data_to_process} | gzip -c >> \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq.gz
            fi
            if [ ! \$(ls ${params.general_folder}*/*/fastq_pass/\${par_basis}/*.fastq | wc -l) -eq 0 ]
            then
            cat \${data_to_process} >> \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq
            fi
        fi      
    }
    export -f fastqmergeParallel 
    start=\$(date +%s)
    if [ -f ${params.run_dir}corruption_cleanup.txt ]
    then
        cleanup=\$(cat ${params.run_dir}corruption_cleanup.txt)
        if [ \$cleanup -eq 1 ]
        then
            echo "0" > ${params.run_dir}corruption_cleanup.txt
        else
            parallel -v -u --env fastqmergeParallel --no-notice -j ${params.threads} fastqmergeParallel ::: ${data_string} ::: ${params.run_dir}
        fi
    else
        parallel -v -u --env fastqmergeParallel --no-notice -j ${params.threads} fastqmergeParallel ::: ${data_string} ::: ${params.run_dir}
    fi
    end=\$(date +%s)
    time=\$(echo "\$((\$end-\$start))")
    echo "fastq_merge,${iteration.value},\$time" >> ${params.run_dir}processing_time_table.csv
    """
    else
    """
    echo "No new data"
    """
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                          Size checkup                                                                   //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process readlength_check{
    input:
    val string from size_checkup
    val string_array from size_checkup2

    script:
    println "Running readlength analysis"
    """
    for i in ${string}
    do
        if [ ${params.barcoded} -eq 1 ]
        then
            basis=\$(basename \$(dirname \${i}))
        else
            basis=\$(basename \$(dirname \$(dirname \$(dirname \${i}))))
            echo \$basis >> ${params.run_dir}error_logs/readlength_counting_error.log
        fi 
        bash ${params.script_dir}get_read_length_from_fastq.sh \$i \$basis ${params.run_dir}ReadLengthFolder || echo ${params.run_dir}error_logs/readlength_counting_error.log
    done 
    """
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                Minimap2 Alignment Genome                                                                //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//Data alignment is performed serial with all new files an iteration of file data_aignment_prep is extracting

process minimap_alignment{
	memory "8 GB"
    
    input:
	val string from minimap_channel
	val string_array from minimap_channel2	
	
    output:
	val bam_string into annotation_channel
	val string_array into annotation_channel2
    
   
	
	script:
	
	
	
	//Alignment running counter restricts number of parallel running alignment series to 2
	bam_string=""
    println "Minimap2 Alignment started"
    sleep(3000)
    for(i in 0..string_array.size()-1){
       if (string_array.size() > 0){
           element = string_array.get(i)
           if (params.barcoded == 1){
           basis = element.getParent().getName()
           }
           else {
               basis = element.getParent().getParent().getParent().getName()
           }
           filename = element.getName()
           converted_filename =  filename.minus(".${params.suffix}")
           outname=params.run_dir + basis + "/bam_files/" + converted_filename + ".bam" + " "
           bam_string=bam_string + outname
       }
    }
    data_string = ""
    for (i in 0..params.sample_names.size()-1){
	    data_string = data_string + params.run_dir + params.sample_names.get(i) + "/bam_files/" + " "
    }
    if (string != "")
	    """
        start=\$(date +%s)
    	for i in ${string}
    	do 
            if [ ${params.barcoded} -eq 1 ]
            then
                basis=\$(basename \$(dirname "\${i}"))
                filename=\$(basename "\${i}")
                converted_filename=\${filename/${params.suffix}/bam}
                outname=${params.run_dir}\${basis}/bam_files/\${converted_filename}
            else
                basis=\$(basename \$(dirname \$(dirname \$(dirname "\${i}"))))
                filename=\$(basename "\${i}")
                converted_filename=\${filename/${params.suffix}/bam}
                outname=${params.run_dir}\${basis}/bam_files/\${converted_filename}
            fi
            if [ ${params.DRS} -eq 1 ]
            then
            minimap2 --MD -ax splice -uf -k14 -t ${params.threads} ${params.run_dir}MT-human_ont.mmi \$i | samtools view -hbS -F 3844 | samtools sort > \${outname} || echo "\$i " >> ${params.run_dir}error_logs/minimap2_genome_failed.log
            else
            minimap2 --MD -ax splice -t ${params.threads} ${params.run_dir}MT-human_ont.mmi \$i | samtools view -hbS -F 3844 | samtools sort > \${outname} || echo "\$i " >> ${params.run_dir}error_logs/minimap2_genome_failed.log 
            fi
        done
        end=\$(date +%s)
        time=\$(echo "\$((\$end-\$start))")
        echo "minimap_genome,${iteration.value},\$time" >> ${params.run_dir}processing_time_table.csv
        function samtoolsParallel  {
        if [ ! \$(ls \${1}| wc -l) -eq 0 ]
        then
            sample=\$1
            run_dir=\$2
            par_basis=\$(basename \$(dirname \${sample}))
            if [ ! -f \${run_dir}bam_genome_merged/\${par_basis}.bam ]
            then
                #samtools merge \${run_dir}\${par_basis}/all_gene.bam \${sample}*.bam -f -h \${sample}*.bam --threads ${params.threads} -c -p || echo "\${sample}" >> \${run_dir}error_logs/merge_genome_all_error.log
                samtools merge \${run_dir}\${par_basis}/all_gene.bam \${sample}*.bam -f --threads ${params.threads} -c -p || echo "\${sample}" >> \${run_dir}error_logs/merge_genome_all_error.log
                cp \${run_dir}\${par_basis}/all_gene.bam \${run_dir}bam_genome_merged/\${par_basis}.bam 
                samtools index \${run_dir}bam_genome_merged/\${par_basis}.bam || echo "Indexing failed" >> \${run_dir}/error_logs/indexing_genome_all_error.log
            else
            bam_files_to_merge=\${run_dir}\${par_basis}/all_gene.bam" "
            for i in ${bam_string}
            do
                if [[ "\$(basename \$(dirname \$(dirname "\${i}")))" == "\${par_basis}" ]]
                then
                    bam_files_to_merge=\${bam_files_to_merge}"\${i}"" "
                fi
            done
            #echo \$bam_files_to_merge >> \${run_dir}/error_logs/bam_files_to_merge.txt
            #samtools merge \${run_dir}bam_genome_merged/\${par_basis}.bam \${bam_files_to_merge} -f -h \${bam_files_to_merge} --threads ${params.threads} -c -p || echo "\${bam_files_to_merge}" >> \${run_dir}error_logs/merge_genome_few_error.log
            samtools merge \${run_dir}bam_genome_merged/\${par_basis}.bam \${bam_files_to_merge} -f --threads ${params.threads} -c -p || echo "\${bam_files_to_merge}" >> \${run_dir}error_logs/merge_genome_few_error.log
            cp \${run_dir}bam_genome_merged/\${par_basis}.bam \${run_dir}\${par_basis}/all_gene.bam
            samtools index \${run_dir}bam_genome_merged/\${par_basis}.bam || echo "Indexing failed" >> \${run_dir}/error_logs/indexing_genome_few_error.log
            fi 
        fi
        }
        export -f samtoolsParallel 
        start=\$(date +%s)
        parallel -v -u --env samtoolsParallel --no-notice -j ${params.threads} samtoolsParallel ::: ${data_string} ::: ${params.run_dir}
        end=\$(date +%s)
        time=\$(echo "\$((\$end-\$start))")
        echo "samtools_genome,${iteration.value},\$time" >> ${params.run_dir}processing_time_table.csv
        echo -e "Sample\tnum_reads\tnum_mapped_reads" > ${params.run_dir}mapping_stats.txt
        for i in \$(ls ${params.run_dir}bam_genome_merged/*.bam)
        do
            echo \$i
            numReads=\$(samtools view -c \$i)
            numMappedReads=\$(samtools view -c -F 260 \$i)
            sample=\$(basename \$i) 
            sample=\${sample/.bam/}
            echo "\$sample\t\$numReads\t\$numMappedReads" >> ${params.run_dir}mapping_stats.txt 
        done
        if [ \$(basename \$(dirname \$(dirname \$(dirname${params.script_dir}))))=="NanopoReaTA" ]
        then
            if [ \$(basename \$(dirname \$(dirname${params.script_dir})))=="app" ]
            then
            rm ${params.script_dir}work/* -r || echo "Erase data in nanopore script directory failed"
            rm ${params.script_dir}.nextflow* -r || echo "Erase nextflow log in script directory failed"
            rm ${params.script_dir}../../work/* -r || echo "Erase data in nanopore work directory failed"
            rm ${params.script_dir}../../.nextflow* -r || echo "Erase nextflow log failed"
            fi
        fi
        """
    else
        """
        echo ""
        """
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                            Minimap2 Alignment Transcriptome                                                             //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


process minimap_transcript_alignment{
    memory "8 GB"
	input:
	val string from minimap_transcript_channel
	val string_array from minimap_transcript_channel2	
	output:
	val bam_string into move_transcripts_channel
	val string_array into move_transcripts_channel2
	
	script:
	sleep(3000)
	
    println "Minimap2 Alignment started"
    bam_string=""
    for(i in 0..string_array.size()-1){
       if (string_array.size() > 0){
           element = string_array.get(i)
           if (params.barcoded == 1){
           basis = element.getParent().getName()
           }
           else {
               basis = element.getParent().getParent().getParent().getName()
           }
           filename=element.getName()
           converted_filename =  filename.minus(".${params.suffix}")
           outname=params.run_dir + basis + "/bam_files_transcripts/" + converted_filename + ".bam" + " "
           bam_string=bam_string + outname
       }
    }

	if (string != "")
        """
        start=\$(date +%s)
    	for i in ${string}
    	do 
            if [ ${params.barcoded} == 1 ]
            then
                basis=\$(basename \$(dirname \${i}))
                filename=\$(basename \${i})
                converted_filename=\${filename/${params.suffix}/bam}
                outname=${params.run_dir}\${basis}/bam_files_transcripts/\${converted_filename}
            else
                basis=\$(basename \$(dirname \$(dirname \$(dirname \${i}))))
                filename=\$(basename \${i})
                converted_filename=\${filename/${params.suffix}/bam}
                outname=${params.run_dir}\${basis}/bam_files_transcripts/\${converted_filename}
            fi
            if [ ${params.DRS} -eq 1 ]
            then
            minimap2 --MD -ax map-ont -uf -k14 -t ${params.threads} ${params.run_dir}MT-human_transcript_ont.mmi ${params.run_dir} \$i | samtools view -hbS -F 3844 | samtools sort > \${outname} || echo "\$i " >> ${params.run_dir}error_logs/minimap2_transcript_failed.log
            else
            minimap2 --MD -ax map-ont -t ${params.threads} ${params.run_dir}MT-human_transcript_ont.mmi ${params.run_dir} \$i | samtools view -hbS -F 3844 | samtools sort > \${outname} || echo "\$i " >> ${params.run_dir}error_logs/minimap2_transcript_failed.log
            fi
        done
        end=\$(date +%s)
        time=\$(echo "\$((\$end-\$start))")
        echo "minimap_transcriptome,${iteration.value},\$time" >> ${params.run_dir}processing_time_table.csv
        """
    else
        """
        echo ""
        """
}

process move_transcript_files{
    memory "4 GB"
    input: 
    val bam_string from move_transcripts_channel
	val string_array from move_transcripts_channel2

    output:
    val bam_full_string into moving_done_channel
    val string_array into moving_done_channel2

    script:
    bam_full_string = ""
    for(i in 0..string_array.size()-1){
       if (string_array.size() > 0){
           element = string_array.get(i)
           if (params.barcoded == 1){
           basis = element.getParent().getName()
           }
           else {
               basis = element.getParent().getParent().getParent().getName()
           }
           filename= element.getName()

           converted_filename =  filename.minus(".${params.suffix}")
           outname= params.run_dir + basis + "/bam_files_transcripts/" + "full/" + converted_filename
           bam_full_string=bam_full_string + outname + ".bam" + " "
       }
    }
    while (salmon_annotation_running.value == 1){
            sleep 4000
    }
    move_transcripts_minimap2_running.value = 1
    if (bam_string != "")
        """
        for i in ${bam_string}
        do
            basis=\$(dirname \${i})
            filename=\$(basename \${i})
            cp \${i} \$basis/full/\$filename || echo "File does not exist: \$basis/full/\$filename" > ${params.run_dir}error_logs/copy_full_bam_failed.log
        done
        """
    else
        """
        echo empty
        """

}

process moving_done{
     input:
     val bam_full_string from moving_done_channel
     val string_array from moving_done_channel2

     output:
     val bam_full_string into salmon_channel
     val string_array into salmon_channel2

     exec:
     if (bam_full_string != ""){
        move_transcripts_minimap2_running.value = 0
        alignment_running_counter_transcript.value = 0
     }
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                 FeatureCount Annotation                                                                 //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


process count_features{
    memory "4 GB"
    input:
    val bam_string from annotation_channel
    val string_array from annotation_channel2
    output:
    val fc_string into merge_all 
    val string_array into merge_all2

    script:
    if (bam_string != ""){
    alignment_running_counter.value = 0
    }
    fc_string = ""
    fc_string_splice = ""
    data_string = ""
    for(i in 0..params.sample_names.size()-1){
       if(string_array.size() > 0){
           basis = params.sample_names.get(i)
           outname= params.run_dir + basis + "/merged_fc/merged_fc.csv" + " "
           fc_string= fc_string + outname
           data_string = data_string + params.sample_names.get(i) + " "
       }
    } 
    if (bam_string != "")
 

    """
    start=\$(date +%s)
    for basis in ${data_string}
    do
        featureCounts -a "${params.genome_gtf}" -F 'GTF' -L -T ${params.threads} -o ${params.run_dir}\${basis}/merged_fc/merged_fc.csv ${params.run_dir}\${basis}/all_gene.bam || echo \${basis} >> ${params.run_dir}error_logs/featureCounts_error.log
    done
    end=\$(date +%s)
    time=\$(echo "\$((\$end-\$start))")
    echo "featureCounts,${iteration.value},\$time" >> ${params.run_dir}processing_time_table.csv
    """
    else
    """
    echo ""
    """
}


//Merge all files from all folders of interest into one feature Count table
process merge_table_of_all_folders{
    memory "4 GB"
    input:
    val fc_string from merge_all
    val string_array from merge_all2
    
    output:
    val fc_string into merging_done_channel
    val string_array into write_processed_data_channel
 

    script:
    println "Feature count merging all files"
    
    //While loop is representing a threadlock situation
    if (fc_string != ""){
        while (fc_merging.value == 1){
            sleep 100
        }
        fc_merging.value = 1
        data_string = ""
        for (i in 0..params.sample_names.size()-1){
	    data_string = data_string + params.run_dir + params.sample_names.get(i) + "/" + " "
        }
        //println data_string
    }
    if (fc_string != "") 
        """
        python ${params.script_dir}merge_all_fc.py ${data_string} ${params.run_dir}merged_all.csv || echo "${fc_string}" >> ${params.run_dir}error_logs/feature_counts_merge_all_error.log
        """
    else
        """
        echo ""
        """

}
  
process fc_merging_done{
    memory "4 GB"
    input:
    val fc_string from merging_done_channel
    val string_array from write_processed_data_channel  

    output:
    val fc_string into fc_done_channel
    val fc_string into start_infer_experiment_processes
    val fc_string into write_data_seen
    val string_array into write_data_seen2
    
    exec:
    //Threadlock enables access again   
        fc_running.value = 0
        fc_merging.value = 0
}


process infer_experiment{
    memory "4 GB"
    input:
    val fc_string from start_infer_experiment_processes

    script:
    if (fc_string != "")
    """
    python ${params.script_dir}infer_experiment_absolute_gene_amount.py -s ${params.run_dir}merged_all.csv -m "" -o ${params.run_dir}exp_genes_counted_per_sample.csv || echo "${fc_string}" >> ${params.run_dir}error_logs/quantify_genes_detected.log
    python ${params.script_dir}infer_experiment_inner_variability.py -s ${params.run_dir}merged_all.csv -m "" -d ${params.run_dir}inner_variability_plot.csv  -o ${params.run_dir}inner_variability_per_sample.csv || echo "${fc_string}" >> ${params.run_dir}error_logs/quantify_inner_variablity_detected.log
    """
    else
    """
    echo empty
    """


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                 Salmon Annotation Transcriptome                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


process salmon_annotation{
    memory "4 GB"
    input:
        val string from salmon_channel
        val string_array from salmon_channel2

    output:
        val string into salmon_annotation_merge
        val string_array into salmon_annotation_merge2

    script:
    if (string != ""){ 
        data_string = ""
        for (i in 0..params.sample_names.size()-1){
	        data_string = data_string + params.run_dir + params.sample_names.get(i) + "/" + "bam_files_transcripts" + "/" + "full" + "/" +  " "
        }
        while (move_transcripts_minimap2_running.value == 1){
            sleep 100
        }
        salmon_annotation_running.value= 1 
    
    bam_string_transcripts=""
    for(i in 0..string_array.size()-1){
       if (string_array.size() > 0){
           element = string_array.get(i)
           if (params.barcoded == 1){
           basis = element.getParent().getName()
           }
           else {
               basis = element.getParent().getParent().getParent().getName()
           }
           filename=element.getName()
           converted_filename =  filename.minus(".${params.suffix}")
           outname=params.run_dir + basis + "/bam_files_transcripts/full/" + converted_filename + ".bam" + " "
           bam_string_transcripts=bam_string_transcripts + outname
       }
    } 
    //println(bam_string_transcripts)
    }

    if (string != "")
        """
        function samtoolsParallel2 {
            if [ ! \$(ls \${1}| wc -l) -eq 0 ] 
            then
                sample=\$1
                run_dir=\$2
                par_basis=\$(basename \$(dirname \$(dirname \${sample})))
                if [ ! -f \${run_dir}bam_transcriptome_merged/\${par_basis}.bam ]
                then
                    #samtools merge \${run_dir}\${par_basis}/salmon/all.bam \${sample}*.bam -f -h \${sample}*.bam --threads ${params.threads} -c -p || echo "\${sample}" >> \${run_dir}/error_logs/merge_transcriptome_all_samtools_error.log
                    samtools merge \${run_dir}\${par_basis}/salmon/all.bam \${sample}*.bam -f --threads ${params.threads} -c -p || echo "\${sample}" >> \${run_dir}/error_logs/merge_transcriptome_all_samtools_error.log
                    cp \${run_dir}\${par_basis}/salmon/all.bam \${run_dir}bam_transcriptome_merged/\${par_basis}.bam 
                    samtools index \${run_dir}bam_transcriptome_merged/\${par_basis}.bam || echo "Indexing failed" >> \${run_dir}/error_logs/indexing_transcriptome_all_error.log
                else
                bam_files_to_merge=\${run_dir}\${par_basis}/salmon/all.bam" "
                for i in ${bam_string_transcripts}
                do
                if [[ "\$(basename \$(dirname \$(dirname \$(dirname \${i}))))" == "\${par_basis}" ]]
                    then
                    bam_files_to_merge=\${bam_files_to_merge}\${i}" "
                    fi
                done
                #echo \$bam_files_to_merge >> \${run_dir}/error_logs/bam_files_to_merge_transcripts.txt
                #samtools merge \${run_dir}bam_transcriptome_merged/\${par_basis}.bam \${bam_files_to_merge} -f -h \${bam_files_to_merge} --threads ${params.threads} -c -p || echo "\${bam_files_to_merge}" >> \${run_dir}/error_logs/merge_transcriptome_few_error.log
                samtools merge \${run_dir}bam_transcriptome_merged/\${par_basis}.bam \${bam_files_to_merge} -f --threads ${params.threads} -c -p || echo "\${bam_files_to_merge}" >> \${run_dir}/error_logs/merge_transcriptome_few_error.log
                cp \${run_dir}bam_transcriptome_merged/\${par_basis}.bam \${run_dir}\${par_basis}/salmon/all.bam
                samtools index \${run_dir}bam_transcriptome_merged/\${par_basis}.bam --threads ${params.threads} || echo "Indexing failed" >> \${run_dir}/error_logs/indexing_transcriptome_few_error.log
                fi
            fi
        }
        export -f samtoolsParallel2 
        start=\$(date +%s)
        parallel -v -u --env samtoolsParallel2 --no-notice -j 5 samtoolsParallel2 ::: ${data_string} ::: ${params.run_dir}
        end=\$(date +%s)
        time=\$(echo "\$((\$end-\$start))")
        echo "samtools_transcriptome,${iteration.value},\$time" >> ${params.run_dir}processing_time_table.csv
        start=\$(date +%s)
        for i in ${data_string}
        do 
            basis=\$(basename \$(dirname \$(dirname \${i})))
            if [ ! \$(ls \${i}| wc -l) -eq 0 ]
            then
            salmon quant -t ${params.transcriptome_fasta} -l A -a ${params.run_dir}\${basis}/salmon/all.bam -o ${params.run_dir}\${basis}/salmon/ -g ${params.genome_gtf} -p ${params.threads} || echo "\$i " >> ${params.run_dir}error_logs/salmon_annotation_fail.log
            fi
        done
        end=\$(date +%s)
        time=\$(echo "\$((\$end-\$start))")
        echo "salmon,${iteration.value},\$time" >> ${params.run_dir}processing_time_table.csv 
        """
    else
        """
        echo ""
        """
}

process merge_salmon_annotation{
    memory "4 GB"
    input:
    val string from salmon_annotation_merge
    val string_array from salmon_annotation_merge2
    
    output:
    val string into salmon_annotation_add_names
    val string_array into salmon_annotation_add_names2

    script:
    println "Salmon merging all files"  
    data_string = ""
    for (i in 0..params.sample_names.size()-1){
	data_string = data_string + params.run_dir + params.sample_names.get(i) + "/" + " "
    }
    //println data_string
    if (string != "") 
        """
        python ${params.script_dir}merge_all_salmon.py ${data_string} ${params.run_dir}salmon_merged_tpm.csv ${params.run_dir}salmon_merged_absolute.csv || echo "${string}" >> ${params.run_dir}error_logs/salmon_merge_error.log
        """
    else
        """
        echo ""
        """
}

process salmon_annotation_done{
    memory "4 GB"
    input:
    val string from salmon_annotation_add_names
    val string_array from salmon_annotation_add_names2

    output:
    val done into salmon_done_channel
    
    exec:
    if (string != ""){
        salmon_annotation_running.value = 0
    }
    done = 1
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                      Reawake next iteration of data analysis                                                            //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


process reawake_next_round{
    input:
    val fc from fc_done_channel
    val salmon from salmon_done_channel
    val fc_string from write_data_seen
    val string_array from write_data_seen2
    val data_string from fastq_merge_done_channel

    exec:
    if (fc_string != ""){
    println "Feature Count merge done for current iteration"
    File file = new File("${params.run_dir}data_seen.txt")
    if (string_array.size() > 0){
        for (i in 0..string_array.size()-1){
            lines = file.readLines()
            file.write("${lines.get(0)}${string_array.get(i)} ")
        }
      }
    }
    if (fc_string == ""){
    sleep(30000)
    }
    def continue_text = new File("${params.run_dir}process_running.txt").getText()
    println "Continue ?"
    println continue_text
    println "____________"
    continue_bool.value = continue_text.toInteger()
    if (continue_bool.value == 1){
        checkup << Channel.from(1)
    }
    else {
        def continue_text2 = new File("${params.run_dir}process_running.txt")
        continue_text2.write "2"
        while(continue_bool.value != 1){
            def continue_text3 = new File("${params.run_dir}process_running.txt").getText()
            continue_bool.value = continue_text3.toInteger()
            if ( continue_bool.value == 1){
                println "Going to continue"
            }
            else {
                sleep(2000)
            }
        }
        checkup << Channel.from(1)
    }
}
