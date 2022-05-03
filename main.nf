#! /usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
            nextflow run main.nf --fasta input.fasta
            nextflow run main.nf --pdb input.pdb


        Mandatory input:
        --fasta                        Sequence fasta file
          OR
        --pdb                          Structure PDB file
        If a sequence fasta file is provided the structure will be predicted using AlphaFold.
        If a pdb file is provided the sequence will be extracted from the structure.

        After intallation the following arguments should be set in nexflow.config before first use:
        --AFPAP_PATH                   Path to AFPAP root directory
        --Pfam_PATH                    Path to Pfam-A.hmm file

        Optional arguments:
        --outputDir                    Output directory (default './output')
        --ligands                      List of ligands in .pdb or .mol2 formats (default false)
                                       ex. --ligands "path/to/ligand1.pdb path/to/ligand2.mol2"
        --md_exhaustiveness            Moleculat docking - numbers of Monte Carlo runs (default 4)
        --md_box_size                  Moleculat docking - grid box x,y,z points (default 20)
        --md_spacing                   Moleculat docking - spacing between grid points (default 0.375)
        --skipAlphaFold                Skip AlphaFold structure prediction (default true)        
        --skipSequenceAnalysis         Skip sequence analysis (default false)
        --skipPfamSearch               Skip Pfam search (default false)
        --skipStructureViewer          Skip raport 3D structure viewer (default false)
        --skipPocketPrediction         Skip pocket prediction (default false)
        --skipPointMutations           Skip point mutations effect prediction (default false)
        --skipMolecularDocking         Skip molecular docking (default false)
        --skipMultiQC                  Skip MultiQC raport generation (default false)
        --help                         This usage statement
        """
}

process configurePipeline {
    output:
        path("AFPAP_output")

    script:
        """
        rm -rf "AFPAP_output/work"
        mkdir -p "AFPAP_output/work/multiqc_files"
        if [ "$params.skipMolecularDocking" != "true" ] && [ "$params.ligands" != "false" ] ; then
            mkdir -p "AFPAP_output/work/docking"
        fi
        """
}

process getFASTAfromPDB {
    input:
        path pdbFile
        path outDir
    output:
        path "$outDir/${pdbFile.baseName}.fasta"
    script:
    """
    python "$params.AFPAP_PATH/bin/AFPAP_FASTA_from_PDB.py" -i $pdbFile -n ${pdbFile.baseName} -o $outDir

    """
}

process sequenceAnalysis {
    input:
        path fastaFile
        path outDir
    output:
        val 0
    script:
    """
    python "$params.AFPAP_PATH/bin/AFPAP_sequence_analysis.py" -i $fastaFile -o $outDir --AFPAPpath $params.AFPAP_PATH
    """
}

process PfamAnalysis {
    input:
        path fastaFile
        path outDir
    output:
        val 0
    script:
    """
    pfam_scan.pl -clan_overlap -align -json pretty -fasta $fastaFile -dir $params.Pfam_PATH > "$outDir/work/pfam.json"
    python "$params.AFPAP_PATH/bin/AFPAP_pfam.py" -j "$outDir/work/pfam.json" -o $outDir --AFPAPpath $params.AFPAP_PATH
    """
}

process prepareStructure {
    input:
        path pdbFile
        path outDir
    output:
        val 0
    script:
    """
    python "$params.AFPAP_PATH/bin/AFPAP_secondary_structure.py" -i $pdbFile -o $outDir --AFPAPpath $params.AFPAP_PATH
    if [ !$params.skipStructureViewer ] ; then
        python "$params.AFPAP_PATH/bin/AFPAP_structure_analysis.py" -i "$outDir/work/proteinStructure.pdb" -o $outDir --AFPAPpath $params.AFPAP_PATH
    fi
    if [ "$params.skipMolecularDocking" != "true" ] && [ "$params.ligands" != "false" ] ; then
        cp "$outDir/work/proteinStructure.pdb" "$outDir/work/docking/receptor.pdb"
        cd "$outDir/work/docking"
        prepare_receptor -r receptor.pdb -o receptor.pdbqt -A "hydrogens"
    fi
    """
}

process PointMutations {
    input:
        path outDir
        val structConfirm
    output:
        val 0
    script:
    """
    python "$params.AFPAP_PATH/bin/AFPAP_point_mutations.py" -i "$outDir/work/proteinStructure.pdb" -o $outDir --AFPAPpath $params.AFPAP_PATH
    """
}

process pocketPrediction {
    input:
        path outDir
        val structConfirm
    output:
        val 0
    script:
    """
    prank predict -f "$outDir/work/proteinStructure.pdb" -o "$outDir/work" -c alphafold
    python "$params.AFPAP_PATH/bin/AFPAP_p2rank_visualization.py" -p "$outDir/work/visualizations/proteinStructure.pdb.pml" -c "$outDir/work/proteinStructure.pdb_predictions.csv" -o $outDir --AFPAPpath $params.AFPAP_PATH
    (cd "$outDir/work/visualizations"; pymol -cq proteinStructure.pdb.pml)
    python "$params.AFPAP_PATH/bin/AFPAP_p2rank_gallery.py" -o $outDir --AFPAPpath $params.AFPAP_PATH 
    """
}

process molecularDocking {
    input:
        path outDir
        each path(ligandFile)
        val structConfirm
        val ppConfirm
    output:
        val 0
    script: 
    if(!ligandFile.exists()){
        log.warn """
        Ligand file does not exist $ligandFile
        """
        """echo 1"""
    }

    else {
        """
        cp $ligandFile "$outDir/work/docking/"
        mkdir -p "$outDir/work/docking/${ligandFile.baseName}"
        prepare_ligand -l $ligandFile -o "$outDir/work/docking/${ligandFile.baseName}/${ligandFile.baseName}.pdbqt"
        python "$params.AFPAP_PATH/bin/AFPAP_molecular_docking.py" -r "$outDir/work/docking/receptor.pdbqt" -l "$outDir/work/docking/${ligandFile.baseName}/${ligandFile.baseName}.pdbqt" -n "${ligandFile.baseName}" -e $params.md_exhaustiveness --box_size $params.md_box_size --spacing $params.md_spacing -o $outDir --AFPAPpath $params.AFPAP_PATH
        cd "$outDir/work/docking/${ligandFile.baseName}"
        pymol -cq generate_complex.pml
        """
    }

}

process MultiQCreport {
    publishDir "${params.outputDir}", mode: 'copy'
    input:
        path outDir
        val saConfirm
        val pfamConfirm
        val structConfirm
        val pocketConfirm
        val mutationConfirm
        val mdConfirm
    output:
        path("./AFPAP_output")
    script:
    """
    echo $saConfirm $pfamConfirm $structConfirm $pocketConfirm $mutationConfirm $mdConfirm >> $outDir/work/ah.txt
    if [ !$params.skipMultiQC ] ; then
        multiqc -f -c "$params.AFPAP_PATH/config/multiqc_config.yaml" --custom-css-file "$params.AFPAP_PATH/config/multiqc_custom_css.css" -o $outDir $outDir    
    fi
    """
}

workflow {
    // Show help message
    if(params.help) {
        helpMessage()
        exit 0
    }

    if(!params.fasta && !params.pdb){
        log.error"Missing input parameters!"
        exit 1
    }
    configurePipeline()
    pipeline_ch = configurePipeline.out

    if(!params.fasta && params.pdb){
        pdb_ch = Channel.fromPath(params.pdb)
        getFASTAfromPDB(pdb_ch, pipeline_ch)
        fasta_ch = getFASTAfromPDB.out
    }
    else if(!params.pdb && params.fasta){
        log.error"No PDB"
        exit 1
    }
    else {
        fasta_ch = Channel.fromPath(params.fasta)
        pdb_ch = Channel.fromPath(params.pdb)
    }
    
    if(!params.skipSequenceAnalysis){
        sequenceAnalysis(fasta_ch, pipeline_ch)
        sa_ch = sequenceAnalysis.out
    }
    else{
        sa_ch = Channel.from(0)
    }
    
    if(!params.skipPfamSearch){
        PfamAnalysis(fasta_ch, pipeline_ch)
        pfam_ch = PfamAnalysis.out
    }
    else{
        pfam_ch = Channel.from(0)
    }
    
    prepareStructure(pdb_ch, pipeline_ch)
    struct_ch = prepareStructure.out

    if(!params.skipPointMutations){
        PointMutations(pipeline_ch, struct_ch)
        mut_ch = PointMutations.out
    }
    else{
        mut_ch = Channel.from(0)
    }

    if(!params.skipPocketPrediction){
        pocketPrediction(pipeline_ch, struct_ch)
        pp_ch = pocketPrediction.out
    }
    else{
        pp_ch = Channel.from(0)
    }

    if(!params.skipMolecularDocking && params.ligands){
        ligands_ch = Channel.from(params.ligands.tokenize()).flatMap{ files(it) }.collect()
        molecularDocking(pipeline_ch, ligands_ch, struct_ch, pp_ch)
        md_ch = molecularDocking.out.collect()
    }
    else{
        md_ch = Channel.from(0)
    }

    MultiQCreport(pipeline_ch, sa_ch, pfam_ch, struct_ch, pp_ch, mut_ch, md_ch)

}