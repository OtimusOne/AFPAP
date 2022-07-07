#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

def helpMessage() {
    log.info """
----------------------------------------------------------------------------------------------------
 AFPAP v1.0
----------------------------------------------------------------------------------------------------
    Usage:
    The typical command for running the pipeline is as follows:
        nextflow run main.nf --fasta input.fasta
        nextflow run main.nf --pdb input.pdb

    Mandatory input:
    --fasta <path>                  Sequence FASTA file
        OR
    --pdb <path>                    Structure PDB file
    If a sequence FASTA file is provided the structure will be predicted using AlphaFold.
    If a PDB file is provided the sequence will be extracted from the structure.

    Optional arguments:
    --outputDir <path>              Output directory (default './output')

    --pfam_path <path>              Path to the Pfam-A.hmm file directory
                                    Should be set in nexflow.config before first use

    --blastdb_path <path>           Path to the BLAST DB used for sequence conservation
                                    Should be set in nexflow.config before first use

    --blastArgs <args>              Arguments sent to PSI-BLAST
                                    (default "-evalue 1e-5 -num_threads 8 -num_iterations 3")

    --colabfoldArgs <args>          Agruments sent to ColabFold
                                    (default "--amber --use-gpu-relax --templates --num-recycle 3
                                              --num-models 5")

    --ligands <path>                List of ligands in .pdb or .mol2 formats (default false)
                                    ex. --ligands "path/to/ligand1.pdb path/to/ligand2.mol2"

    --dock_pockets <0/1>            Molecular docking - dock ligand against predicted pockets,
                                                        if available (default true)

    --md_exhaustiveness <int>       Molecular docking - numbers of Monte Carlo runs (default 4)

    --md_box_size <int>             Molecular docking - nr. of grid box x,y,z points (default 20)

    --md_spacing <float>            Molecular docking - spacing between grid points (default 0.375)

    --pdb_type <0/1/2>              PDB file type, AlphaFold - 0, X-ray - 1, other - 2 (default 0)

    --skipSequenceProperties <0/1>  Skip sequence properties (default false)

    --skipPfamSearch <0/1>          Skip Pfam search (default false)

    --skipConservationMSA <0/1>     Skip Conservation & Multiple Sequence Alignment (default false)

    --skipStructuralAnalysis <0/1>  Skip structural analysis (default false)

    --skipAlphaFold <0/1>           Skip AlphaFold structure prediction (default false)

    --skipStructureViewer <0/1>     Skip report iCn3D structure viewer (default false)

    --skipPocketPrediction <0/1>    Skip pocket prediction (default false)

    --skipStabilityChanges <0/1>    Skip residue substitution effect prediction (default false)

    --skipMolecularDocking <0/1>    Skip molecular docking (default false)

    --skipMultiQC <0/1>             Skip MultiQC report generation (default false)

    --help                          This usage statement
    """
}

def validateParameters() {
    if (!params.fasta && !params.pdb) {
        log.error'Missing input file parameters!'
        exit 1
    }
    if (!(params.md_exhaustiveness instanceof java.lang.Integer)) {
        log.error'Exhaustiveness must be integer!'
        exit 1
    }
    if (!(params.md_box_size instanceof java.lang.Integer)) {
        log.error'Box size must be integer!'
        exit 1
    }
    if (!(params.pdb_type instanceof java.lang.Integer)) {
        log.error'pdb_type must be integer!'
        exit 1
    }
    if (!params.md_spacing.toString().isNumber()) {
        log.error'Grid spacing must be a number!'
        exit 1
    }

    validBool = [0, 1, false, true, 'false', 'true']
    if (!validBool.contains(params.dock_pockets)) {
        log.error'dock_pockets must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipAlphaFold)) {
        log.error'skipAlphaFold must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipSequenceProperties)) {
        log.error'skipSequenceProperties must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipPfamSearch)) {
        log.error'skipPfamSearch must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipConservationMSA)) {
        log.error'skipConservationMSA must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipStructuralAnalysis)) {
        log.error'skipStructuralAnalysis must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipStructureViewer)) {
        log.error'skipStructureViewer must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipStabilityChanges)) {
        log.error'skipStabilityChanges must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipPocketPrediction)) {
        log.error'skipPocketPrediction must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipMolecularDocking)) {
        log.error'skipMolecularDocking must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipMultiQC)) {
        log.error'skipMultiQC must be 0/1!'
        exit 1
    }
    paramList = ['pfam_path', 'blastdb_path', 'fasta', 'pdb', 'outputDir', 'output-dir', 'help', 'ligands', 'dock_pockets', 'md_exhaustiveness', 'md_box_size', 'md_spacing', 'blastArgs', 'blast-args', 'colabfoldArgs', 'colabfold-args', 'pdb_type', 'skipAlphaFold', 'skip-alpha-fold', 'skipStructuralAnalysis', 'skip-structural-analysis', 'skipSequenceProperties', 'skip-sequence-properties', 'skipPfamSearch', 'skip-pfam-search', 'skipConservationMSA', 'skip-conservation-MSA', 'skipStructureViewer', 'skip-structure-viewer', 'skipStabilityChanges', 'skip-stability-changes', 'skipMolecularDocking', 'skip-molecular-docking', 'skipMultiQC', 'skip-multi-QC', 'skip-pocket-prediction', 'skipPocketPrediction']
    for (parameter in params) {
        if (!paramList.contains(parameter.key)) {
            log.warn"Unknown parameter ${parameter.key}..."
        }
    }
}

// Configure pipeline directory structure
process configurePipeline {
    input:
        path inputFile
    output:
        path("./$params.inputBaseName")

    script:
    params.inputBaseName = "${inputFile.baseName}"
        """
        rm -rf "$params.inputBaseName/work"
        mkdir -p "$params.inputBaseName/work/multiqc_files"
        if ( [ "$params.skipMolecularDocking" != "true" ] && [ "$params.skipMolecularDocking" != 1 ] ) && ( [ "$params.ligands" != "false" ] && [ "$params.ligands" != 0 ] ) ; then
            mkdir -p "$params.inputBaseName/work/docking"
        fi
        if ( [ "$params.skipConservationMSA" != "true" ] && [ "$params.skipConservationMSA" != 1 ] ) ; then
            mkdir -p "$params.inputBaseName/work/MSA"
        fi
        cp $inputFile "./$params.inputBaseName/"
        """
}

// Extract and validate protein sequence from input file
process validateSequence {
    input:
        path pdbFile
        path outDir
        val fileType
    output:
        path "$outDir/${params.inputBaseName}.fasta"
    script:
    """
    python "$projectDir/bin/AFPAP_validate_sequence.py" -i $pdbFile -n $params.inputBaseName -o $outDir -t $fileType
    """
}

// Predict PDB structure using ColabFold
process predictPDB {
    input:
        path fastaFile
        path outDir
    output:
        path "$outDir/${params.inputBaseName}.pdb"
    script:
    """
    colabfold_batch $params.colabfoldArgs $fastaFile $outDir/work/alphafold

    if [[ -n `find $outDir/work/alphafold -name "*_relaxed_rank_1*.pdb" -print -quit` ]]
    then
        cp `find $outDir/work/alphafold -name "*_relaxed_rank_1*.pdb" -print -quit` "$outDir/${params.inputBaseName}.pdb"
    else
        cp `find $outDir/work/alphafold -name "*_unrelaxed_rank_1*.pdb" -print -quit` "$outDir/${params.inputBaseName}.pdb"
    fi
    """
}

// Generate sequence viewer and properties
process sequenceProperties {
    input:
        path fastaFile
        path outDir
    output:
        val 0
    script:
    """
    python "$projectDir/bin/AFPAP_sequence_properties.py" -i $fastaFile -o $outDir --AFPAPpath $projectDir
    """
}

// Align sequence using the Pfam database
process pfamSearch {
    input:
        path fastaFile
        path outDir
    output:
        val 0
    script:
    """
    pfam_scan.pl -clan_overlap -align -json pretty -fasta $fastaFile -dir $params.pfam_path > "$outDir/work/${params.inputBaseName}_pfam.json"
    python "$projectDir/bin/AFPAP_pfam.py" -j "$outDir/work/${params.inputBaseName}_pfam.json" -o $outDir --AFPAPpath $projectDir
    """
}

// Generate MSA and calculate sequence conservation
process conservationMSA {
    input:
        path fastaFile
        path outDir
    output:
        val 0
    script:
    """
    blastHits="$outDir/work/MSA/${params.inputBaseName}.hits"
    blastSequences="\$(tempfile)"
    blastResults="$outDir/work/MSA/${params.inputBaseName}.blast"
    modifiedInputFile="\$(tempfile)"
    blastTempFile="\$(tempfile)"
    muscleTempFile="\$(tempfile)"
    muscleAlignment="$outDir/work/MSA/${params.inputBaseName}.muscle"
    conservationFile="$outDir/work/MSA/${params.inputBaseName}.conservation"

    # Run PSI-BLAST to find ids all similar sequences.
    psiblast < $fastaFile -db $params.blastdb_path $params.blastArgs -outfmt '6 sallseqid qcovs pident' | $projectDir/bin/MSA_filter_BLAST.awk > "\${blastHits}"

    # Get full sequences.
    blastdbcmd -db $params.blastdb_path -entry_batch "\${blastHits}" > "\${blastSequences}"

    # Filter using CD-HIT.
    cd-hit -i "\${blastSequences}" -o "\${blastResults}" >&2

    numSeq=`grep < "\${blastResults}" '^>' | wc -l`
    echo Found "\${numSeq}" sequences in $params.blastdb_path

    if (( "\${numSeq}" >= 50 )); then
        # Change the description from the file to find it later.
        sed < $fastaFile 's/^>/>input_sequence|/' > "\${modifiedInputFile}"

        # Run muscle. Note we need to concat the query sequence in order to get its conservation later.
        cat "\${modifiedInputFile}" "\${blastResults}" > "\${blastTempFile}"
        muscle -align "\${blastTempFile}" -output "\${muscleTempFile}"
        awk -f $projectDir/bin/MSA_sort_MUSCLE.awk < "\${muscleTempFile}" > "\${muscleAlignment}"

        python $projectDir/bin/score_conservation.py -m $projectDir/config/blosum62.bla "\${muscleAlignment}" > "\${conservationFile}"
        cp "\${conservationFile}" "$outDir/work/multiqc_files/"
    else
        echo "Too few hits, skipping alignment..."
    fi
    rm "\${modifiedInputFile}"
    rm "\${blastSequences}"
    rm "\${blastTempFile}"
    rm "\${muscleTempFile}"
    """
}

// Add secondary structure, prepare protein for molecular docking
process prepareStructure {
    input:
        path pdbFile
        path outDir
    output:
        val 0
    script:
    """
    python "$projectDir/bin/AFPAP_secondary_structure.py" -i "$outDir/${params.inputBaseName}.pdb" -o $outDir --name "${params.inputBaseName}_fixed.pdb" --AFPAPpath $projectDir
    if [ "$params.skipStructureViewer" != "true" ] && [ "$params.skipStructureViewer" != 1 ] ; then
        python "$projectDir/bin/AFPAP_structure_viewer.py" -i "$outDir/work/${params.inputBaseName}_fixed.pdb" -o $outDir --AFPAPpath $projectDir --pdb_type $params.pdb_type
    fi
    if ( [ "$params.skipMolecularDocking" != "true" ] && [ "$params.skipMolecularDocking" != 1 ] ) && ( [ "$params.ligands" != "false" ] && [ "$params.ligands" != 0 ] ) ; then
        cp "$outDir/work/${params.inputBaseName}_fixed.pdb" "$outDir/work/docking/receptor.pdb"
        cd "$outDir/work/docking"
        prepare_receptor -r receptor.pdb -o receptor_temp.pdbqt -e
        prepare_receptor -r receptor_temp.pdbqt -o receptor.pdbqt -A "hydrogens" -e
        rm receptor_temp.pdbqt receptor.pdb
    fi
    """
}

// Predict amino acid substitution effect on protein stability
process stabilityChanges {
    input:
        path outDir
        val structConfirm
    output:
        val 0
    script:
    """
    python "$projectDir/bin/AFPAP_stability_changes.py" -i "$outDir/work/${params.inputBaseName}_fixed.pdb" -o $outDir --AFPAPpath $projectDir
    """
}

// Predict binding pockets using P2Rank
process pocketPrediction {
    input:
        path outDir
        val structConfirm
    output:
        val 0
    script:
    if (params.pdb_type != 1) {
        predictionMode = '-c alphafold'
    }
    else {
        predictionMode = ''
    }
    """
    prank predict -f "$outDir/work/${params.inputBaseName}_fixed.pdb" -o "$outDir/work" $predictionMode
    python "$projectDir/bin/AFPAP_p2rank_visualization.py" -p "$outDir/work/visualizations/${params.inputBaseName}_fixed.pdb.pml" --csv "$outDir/work/${params.inputBaseName}_fixed.pdb_predictions.csv" -o $outDir --AFPAPpath $projectDir
    (cd "$outDir/work/visualizations"; pymol -cq "${params.inputBaseName}_fixed.pdb.pml")
    python "$projectDir/bin/AFPAP_p2rank_gallery.py" -o $outDir --AFPAPpath $projectDir
    mv "$outDir/work/run.log" "$outDir/work/visualizations/run.log"
    mv "$outDir/work/params.txt" "$outDir/work/visualizations/params.txt"
    """
}

// Dock ligands using AutoDock Vina
process molecularDocking {
    input:
        path outDir
        each path(ligandFile)
        val structConfirm
        val ppConfirm
    output:
        val 0
    script:
    if (ligandFile.exists()) {
        if (params.dock_pockets) {
            dockMode = '--dock_pockets'
        }
        else {
            dockMode = ''
        }
        """
        cp $ligandFile "$outDir/work/docking/"
        mkdir -p "$outDir/work/docking/${ligandFile.baseName}"
        prepare_ligand -l $ligandFile -o "$outDir/work/docking/${ligandFile.baseName}/${ligandFile.baseName}.pdbqt"
        python "$projectDir/bin/AFPAP_molecular_docking.py" -r "$outDir/work/docking/receptor.pdbqt" -l "$outDir/work/docking/${ligandFile.baseName}/${ligandFile.baseName}.pdbqt" -n "${ligandFile.baseName}" $dockMode --exhaustiveness $params.md_exhaustiveness --box_size $params.md_box_size --spacing $params.md_spacing -o $outDir --AFPAPpath $projectDir
        cd "$outDir/work/docking/${ligandFile.baseName}"
        pymol -cq generate_complex.pml
        """
    }
    else {
        log.warn """
        Ligand file does not exist: $ligandFile
        """
        """echo Ligand file does not exist: $ligandFile"""
    }
}

// Generate MultiQC report and publish results
process processPipelineOutput {
    publishDir "${params.outputDir}", mode: 'copy'
    input:
        path outDir
        val saConfirm
        val pfamConfirm
        val msaConfirm
        val structConfirm
        val pocketConfirm
        val mutationConfirm
        val mdConfirm
    output:
        path("./$params.inputBaseName")
    script:
    """
    if [ "$params.skipMultiQC" != "true" ] && [ "$params.skipMultiQC" != 1 ] ; then
        multiqc -f -c "$projectDir/config/multiqc_config.yaml" --custom-css-file "$projectDir/config/multiqc_custom_css.css" -o $outDir "$outDir/work/multiqc_files"
    fi
    """
}

workflow {
    // --------- Pipeline Configuration ---------

    // Show help message
    if (params.help) {
        helpMessage()
        exit 0
    }

    validateParameters()

    if (params.pdb) {
        pdb_ch = Channel.fromPath(params.pdb)
        configurePipeline(pdb_ch)
        pipeline_ch = configurePipeline.out
        validateSequence(pdb_ch, pipeline_ch, 1)
        fasta_ch = validateSequence.out
    }
    else {
        raw_fasta_ch = Channel.fromPath(params.fasta)
        configurePipeline(raw_fasta_ch)
        pipeline_ch = configurePipeline.out
        validateSequence(raw_fasta_ch, pipeline_ch, 0)
        fasta_ch = validateSequence.out

        if (params.skipStructuralAnalysis) {
            pdb_ch = Channel.from(0)
        }
        else if (params.skipAlphaFold) {
            log.error'No AlphaFold - provide PDB file...'
            exit 1
        }
        else {
            predictPDB(fasta_ch, pipeline_ch)
            pdb_ch = predictPDB.out
        }
    }

    // --------- Sequence Analysis ---------

    if (!params.skipSequenceProperties) {
        sequenceProperties(fasta_ch, pipeline_ch)
        sa_ch = sequenceProperties.out
    }
    else {
        sa_ch = Channel.from(0)
    }

    if (!params.skipPfamSearch) {
        pfamSearch(fasta_ch, pipeline_ch)
        pfam_ch = pfamSearch.out
    }
    else {
        pfam_ch = Channel.from(0)
    }

    if (!params.skipConservationMSA) {
        conservationMSA(fasta_ch, pipeline_ch)
        msa_ch = conservationMSA.out
    }
    else {
        msa_ch = Channel.from(0)
    }

    if (!params.skipStructuralAnalysis) {
        prepareStructure(pdb_ch, pipeline_ch)
        struct_ch = prepareStructure.out
    }
    else {
        struct_ch = Channel.from(0)
    }

    // --------- Structural Analysis ---------

    if (!params.skipStabilityChanges && !params.skipStructuralAnalysis) {
        stabilityChanges(pipeline_ch, struct_ch)
        mut_ch = stabilityChanges.out
    }
    else {
        mut_ch = Channel.from(0)
    }

    if (!params.skipPocketPrediction && !params.skipStructuralAnalysis) {
        pocketPrediction(pipeline_ch, struct_ch)
        pp_ch = pocketPrediction.out
    }
    else {
        pp_ch = Channel.from(0)
    }

    if (!params.skipMolecularDocking && params.ligands && !params.skipStructuralAnalysis) {
        ligands_ch = Channel.from(params.ligands.tokenize()).flatMap { files(it) }.collect()
        molecularDocking(pipeline_ch, ligands_ch, struct_ch, pp_ch)
        md_ch = molecularDocking.out.collect()
    }
    else {
        md_ch = Channel.from(0)
    }

    // --------- Generate MultiQC report ---------
    processPipelineOutput(pipeline_ch, sa_ch, pfam_ch, msa_ch, struct_ch, pp_ch, mut_ch, md_ch)
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
