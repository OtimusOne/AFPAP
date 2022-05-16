#! /usr/bin/env nextflow
/* groovylint-disable LineLength */
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
    --fasta <path>                  Sequence fasta file
        OR
    --pdb <path>                    Structure PDB file
    If a sequence fasta file is provided the structure will be predicted using AlphaFold.
    If a pdb file is provided the sequence will be extracted from the structure.

    --AFPAP_PATH <path>             Path to AFPAP root directory
                                    Should be set in nexflow.config before first use

    Optional arguments:
    --outputDir <path>              Output directory (default './output')

    --Pfam_PATH <path>              Path to Pfam-A.hmm file.
                                    Should be set in nexflow.config before first use

    --ligands <path>                List of ligands in .pdb or .mol2 formats (default false)
                                    ex. --ligands "path/to/ligand1.pdb path/to/ligand2.mol2"

    --dock_pockets <0/1>            Molecular docking - dock ligand against predicted pockets,
                                                        if available (default true)

    --md_exhaustiveness <int>       Molecular docking - numbers of Monte Carlo runs (default 4)

    --md_box_size <int>             Molecular docking - nr. of grid box x,y,z points (default 20)

    --md_spacing <float>            Molecular docking - spacing between grid points (default 0.375)

    --pdb_type <0/1/2>              PDB file type, AlphaFold - 0, X-ray - 1, other - 2 (default 0)

    --skipSequenceAnalysis <0/1>    Skip sequence analysis (default false)

    --skipPfamSearch <0/1>          Skip Pfam search (default false)

    --skipStructuralAnalysis <0/1>  Skip structual analysis (default false)

    --skipAlphaFold <0/1>           Skip AlphaFold structure prediction (default false)

    --skipStructureViewer <0/1>     Skip raport 3D structure viewer (default false)

    --skipPocketPrediction <0/1>    Skip pocket prediction (default false)

    --skipPointMutations <0/1>      Skip point mutations effect prediction (default false)

    --skipMolecularDocking <0/1>    Skip molecular docking (default false)

    --skipMultiQC <0/1>             Skip MultiQC raport generation (default false)

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
    if (!validBool.contains(params.skipStructuralAnalysis)) {
        log.error'skipStructuralAnalysis must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipAlphaFold)) {
        log.error'skipAlphaFold must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipSequenceAnalysis)) {
        log.error'skipSequenceAnalysis must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipPfamSearch)) {
        log.error'skipPfamSearch must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipStructureViewer)) {
        log.error'skipStructureViewer must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipPocketPrediction)) {
        log.error'skipPocketPrediction must be 0/1!'
        exit 1
    }
    if (!validBool.contains(params.skipPointMutations)) {
        log.error'skipPointMutations must be 0/1!'
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
    paramList = ['AFPAP_PATH', 'Pfam_PATH', 'fasta', 'pdb', 'outputDir', 'output-dir', 'help', 'ligands', 'dock_pockets', 'md_exhaustiveness', 'md_box_size', 'md_spacing', 'colabfoldArgs', 'colabfold-args', 'pdb_type', 'skipAlphaFold', 'skip-alpha-fold', 'skipStructuralAnalysis', 'skip-structural-analysis', 'skipSequenceAnalysis', 'skip-sequence-analysis', 'skipPfamSearch', 'skip-pfam-search', 'skipStructureViewer', 'skip-structure-viewer', 'skipPointMutations', 'skip-point-mutations', 'skipMolecularDocking', 'skip-molecular-docking', 'skipMultiQC', 'skip-multi-QC', 'skip-pocket-prediction', 'skipPocketPrediction']
    for (parameter in params) {
        if (!paramList.contains(parameter.key)) {
            log.warn"Unknown parameter ${parameter.key}..."
        }
    }
}

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
        cp $inputFile "./$params.inputBaseName/"
        """
}

process getFASTAfromPDB {
    input:
        path pdbFile
        path outDir
    output:
        path "$outDir/${params.inputBaseName}.fasta"
    script:
    """
    python "$params.AFPAP_PATH/bin/AFPAP_FASTA_from_PDB.py" -i $pdbFile -n "${params.inputBaseName}.fasta" -o $outDir
    """
}

process predictPDB {
    input:
        path fastaFile
        path outDir
    output:
        path "$outDir/${params.inputBaseName}.pdb"
    script:
    """
    colabfold_batch $params.colabfoldArgs $fastaFile $outDir/work/alphafold

    if [[ -n `find $outDir/work/alphafold -name "*_relaxed_rank_1*.pdb"` ]]
    then
        cp `find $outDir/work/alphafold -name "*_relaxed_rank_1*.pdb"` "$outDir/${params.inputBaseName}.pdb"
    else
        cp `find $outDir/work/alphafold -name "*_unrelaxed_rank_1*.pdb"` "$outDir/${params.inputBaseName}.pdb"
    fi
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

process pfamSearch {
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
    python "$params.AFPAP_PATH/bin/AFPAP_secondary_structure.py" -i "$outDir/${params.inputBaseName}.pdb" -o $outDir --name "${params.inputBaseName}_fixed.pdb" --AFPAPpath $params.AFPAP_PATH
    if [ "$params.skipStructureViewer" != "true" ] && [ "$params.skipStructureViewer" != 1 ] ; then
        python "$params.AFPAP_PATH/bin/AFPAP_structure_viewer.py" -i "$outDir/work/"${params.inputBaseName}_fixed.pdb"" -o $outDir --AFPAPpath $params.AFPAP_PATH --pdb_type $params.pdb_type
    fi
    if ( [ "$params.skipMolecularDocking" != "true" ] && [ "$params.skipMolecularDocking" != 1 ] ) && ( [ "$params.ligands" != "false" ] && [ "$params.ligands" != 0 ] ) ; then
        cp "$outDir/work/${params.inputBaseName}_fixed.pdb" "$outDir/work/docking/receptor.pdb"
        cd "$outDir/work/docking"
        prepare_receptor -r receptor.pdb -o receptor_temp.pdbqt -e
        prepare_receptor -r receptor_temp.pdbqt -o receptor.pdbqt -A "hydrogens" -e
    fi
    """
}

process pointMutations {
    input:
        path outDir
        val structConfirm
    output:
        val 0
    script:
    """
    python "$params.AFPAP_PATH/bin/AFPAP_point_mutations.py" -i "$outDir/work/${params.inputBaseName}_fixed.pdb" -o $outDir --AFPAPpath $params.AFPAP_PATH
    """
}

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
    python "$params.AFPAP_PATH/bin/AFPAP_p2rank_visualization.py" -p "$outDir/work/visualizations/${params.inputBaseName}_fixed.pdb.pml" --csv "$outDir/work/${params.inputBaseName}_fixed.pdb_predictions.csv" -o $outDir --AFPAPpath $params.AFPAP_PATH
    (cd "$outDir/work/visualizations"; pymol -cq "${params.inputBaseName}_fixed.pdb.pml")
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
        python "$params.AFPAP_PATH/bin/AFPAP_molecular_docking.py" -r "$outDir/work/docking/receptor.pdbqt" -l "$outDir/work/docking/${ligandFile.baseName}/${ligandFile.baseName}.pdbqt" -n "${ligandFile.baseName}" $dockMode --exhaustiveness $params.md_exhaustiveness --box_size $params.md_box_size --spacing $params.md_spacing -o $outDir --AFPAPpath $params.AFPAP_PATH
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

process processPipelineOutput {
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
        path("./$params.inputBaseName")
    script:
    """
    echo $saConfirm $pfamConfirm $structConfirm $pocketConfirm $mutationConfirm $mdConfirm >> $outDir/work/ah.txt
    if [ "$params.skipMultiQC" != "true" ] && [ "$params.skipMultiQC" != 1 ] ; then
        multiqc -f -c "$params.AFPAP_PATH/config/multiqc_config.yaml" --custom-css-file "$params.AFPAP_PATH/config/multiqc_custom_css.css" -o $outDir "$outDir/work/" 
    fi
    """
}

workflow {
    // Show help message
    if (params.help) {
        helpMessage()
        exit 0
    }

    validateParameters()

    if (!params.fasta && params.pdb) {
        pdb_ch = Channel.fromPath(params.pdb)
        configurePipeline(pdb_ch)
        pipeline_ch = configurePipeline.out
        getFASTAfromPDB(pdb_ch, pipeline_ch)
        fasta_ch = getFASTAfromPDB.out
    }
    else if (!params.pdb && params.fasta) {
        if (params.skipStructuralAnalysis) {
            fasta_ch = Channel.fromPath(params.fasta)
            configurePipeline(fasta_ch)
            pipeline_ch = configurePipeline.out
            pdb_ch = Channel.from(0)
        }
        else if (params.skipAlphaFold) {
            log.error'No AlphaFold - provide PDB file..'
            exit 1
        }
        else {
            fasta_ch = Channel.fromPath(params.fasta)
            configurePipeline(fasta_ch)
            pipeline_ch = configurePipeline.out
            predictPDB(fasta_ch, pipeline_ch)
            pdb_ch = predictPDB.out
        }
    }
    else {
        fasta_ch = Channel.fromPath(params.fasta)
        pdb_ch = Channel.fromPath(params.pdb)
        configurePipeline(pdb_ch)
        pipeline_ch = configurePipeline.out
    }

    if (!params.skipSequenceAnalysis) {
        sequenceAnalysis(fasta_ch, pipeline_ch)
        sa_ch = sequenceAnalysis.out
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

    if (!params.skipStructuralAnalysis) {
        prepareStructure(pdb_ch, pipeline_ch)
        struct_ch = prepareStructure.out
    }
    else {
        struct_ch = Channel.from(0)
    }

    if (!params.skipPointMutations && !params.skipStructuralAnalysis) {
        pointMutations(pipeline_ch, struct_ch)
        mut_ch = pointMutations.out
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

    processPipelineOutput(pipeline_ch, sa_ch, pfam_ch, struct_ch, pp_ch, mut_ch, md_ch)
}
