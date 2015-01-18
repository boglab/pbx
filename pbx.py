import datetime
import subprocess
import os
import os.path
import ConfigParser
import shlex
import sys

def run_command(options, command):
    #subprocess.check_call(["bash", "-x", "-e"] + shlex.split(command))
    subprocess.check_call(["bash", "-e"] + shlex.split(command))

def validate_options(options):
    pass

def get_utils_and_citations():
    return """
    The PBX toolkit is built on free software. Please cite the following papers for this run:
    
    AMOS:
    Treangen TJ, Sommer DD, Angly FE, Koren S, Pop M: Next Generation Sequence Assembly with AMOS. In Current Protocols in Bioinformatics. John Wiley & Sons, Inc.; 2011
    
    BLASR:
    Chaisson MJ, Tesler G: Mapping single molecule sequencing reads using basic local alignment with successive refinement (BLASR): application and theory. Bmc Bioinformatics 2012, 13.
    
    BLAST+:
    Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden T: BLAST+: architecture and applications. BMC Bioinformatics 2009, 10:421.
    
    Celera Assembler:
    Myers EW, Sutton GG, Delcher AL, Dew IM, Fasulo DP, Flanigan MJ, Kravitz SA, Mobarry CM, Reinert KHJ, Remington KA, et al: A whole-genome assembly of Drosophila. Science 2000, 287:2196-2204.
    
    HGAP:
    Chin CS, Alexander DH, Marks P, Klammer AA, Drake J, Heiner C, Clum A, Copeland A, Huddleston J, Eichler EE, et al: Nonhybrid, finished microbial genome assemblies from long-read SMRT sequencing data. Nature Methods 2013, 10:563-+.
    
    Minimus:
    Sommer DD, Delcher AL, Salzberg SL, Pop M. Minimus: a fast, lightweight genome assembler. BMC Bioinformatics. 2007 Feb 26;8:64.
    
    MUMmer:
    Kurtz S, Phillippy A, Delcher AL, Smoot M, Shumway M, Antonescu C, Salzberg SL: Versatile and open software for comparing large genomes. Genome Biology 2004, 5.
    """

def get_options(base_path):
    
    config = ConfigParser.RawConfigParser()
    config.read(sys.argv[1])
    
    options = {}
    
    options["stop_after"] = "export"
    
    if config.has_option('General', 'stop_after'):
        options["stop_after"] = config.get('General', 'stop_after')
    
    options["dry_run"] = False
    
    if config.has_option('General', 'dry_run'):
        options["dry_run"] = config.getboolean('General', 'dry_run')
    
    # Environment
    
    options["env"] = {
        # This is needed to stop SMRTAnalysis 2.3 from removing our environment variables when we source setup.sh
        "SMRT_ENV_PASSTHROUGH_VARS": "PBX_.*",
        # Auto
        "PBX_BASE_PATH": base_path,
        "PBX_SCRIPTS_PATH": os.path.join(base_path, "scripts"),
        "PBX_PROTOCOL_TEMPLATES_PATH": os.path.join(base_path, "protocol_templates"),
        "PBX_WORKFLOW_PATH": os.path.join(base_path, "workflow"),
        # User config
        # General
        "PBX_SMRTANALYSIS_PATH": config.get('General', 'smrtanalysis_path'),
        "PBX_MUMMER_PATH": config.get('General', 'mummer_path'),
        "PBX_RAW_READS_PATH": config.get('General', 'raw_reads_path'),
        # Whitelisting
        "PBX_TALE_SEQS_WHITELISTING": os.path.join(base_path, "tale_seqs", "whitelisting", config.get('Whitelisting', 'tale_seqs_file_whitelisting')),
        # Preassembly
        "PBX_PREASSEMBLY_MIN_SEED_READ_LENGTH": config.get('Preassembly', 'min_seed_read_length'),
        "PBX_PREASSEMBLY_MIN_SUBREAD_LENGTH": config.get('Preassembly', 'min_subread_length'),
        "PBX_PREASSEMBLY_MIN_TRIMMED_PREASSEMBLED_READ_QV": config.get('Preassembly', 'min_trimmed_preassembled_read_qv'),
        "PBX_PREASSEMBLY_MIN_TRIMMED_PREASSEMBLED_READ_LENGTH": config.get('Preassembly', 'min_trimmed_preassembled_read_length'),
        # Resequencing
        "PBX_RESEQUENCING_MIN_SUBREAD_LENGTH": config.get('Resequencing', 'min_subread_length'),
        # Export
        "PBX_TALE_SEQS_EXPORT": os.path.join(base_path, "tale_seqs", "exporter", config.get('Export', 'tale_seqs_file_export')),
        "PBX_TALE_SEQS_EXPORT_BOUNDARIES": os.path.join(base_path, "tale_seqs", "exporter", config.get('Export', 'tale_seqs_file_export_boundaries')),
    }
    
    if config.has_option('General', 'results_path'):
        options["env"]["PBX_RESULTS_PATH"] = config.get('General', 'results_path')
    else:
        options["env"]["PBX_RESULTS_PATH"] = os.path.join(base_path, "results", "results_" + datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S'))

    # Other auto set options
    
    options["env"]["PBX_WHITELISTING_RESULTS_PATH"] = os.path.join(options["env"]["PBX_RESULTS_PATH"], "blasr_raw_to_bls_tals_alignments")
    options["env"]["PBX_PREASSEMBLY_RESULTS_PATH"] = os.path.join(options["env"]["PBX_RESULTS_PATH"], "preassembly_job")
    options["env"]["PBX_ASSEMBLY_RESULTS_PATH"] = os.path.join(options["env"]["PBX_RESULTS_PATH"], "minimo_results")
    options["env"]["PBX_RESEQUENCING_RESULTS_PATH"] = os.path.join(options["env"]["PBX_RESULTS_PATH"], "resequencing")
    options["env"]["PBX_RESEQUENCED_CONTIG_ASSEMBLY_RESULTS_PATH"] = os.path.join(options["env"]["PBX_RESULTS_PATH"], "combine_resequenced_tals")
    
    return options

if __name__ == "__main__":
    
    base_path = os.path.dirname(os.path.realpath(__file__))
    
    options = get_options(base_path)
    
    os.environ.update(options["env"])
    
    print(get_utils_and_citations())
    
    print("Identifying TALE region reads...")
    
    run_command(options, "%s/run_tale_whitelisting.sh" % options["env"]["PBX_WORKFLOW_PATH"])
    
    if options["stop_after"] == "whitelisting":
        exit()
    
    print("Correcting errors in long TALE region reads using short reads...")
    
    run_command(options, "%s/run_tale_preassembly.sh" % options["env"]["PBX_WORKFLOW_PATH"])
    
    if options["stop_after"] == "preassembly":
        exit()
    
    print("Assembling error-corrected TALE region reads with various settings...")
    
    run_command(options, "%s/run_tale_read_assembly.sh" % options["env"]["PBX_WORKFLOW_PATH"])
    
    if options["stop_after"] == "assembly":
        exit()
    
    print("Polishing assembled TALE regions... (this will take several hours. Check back tomorrow.")
    
    run_command(options, "%s/run_tale_resequencing.sh" % options["env"]["PBX_WORKFLOW_PATH"])
    
    if options["stop_after"] == "resequencing":
        exit()
    
    print("Merging polished TALE regions.")
    
    run_command(options, "%s/run_tale_resequenced_contig_assembly.sh" % options["env"]["PBX_WORKFLOW_PATH"])
    
    if options["stop_after"] == "resequenced_contig_assembly":
        exit()
    
    print("Exporting TALE RVD sequences from polished assemblies. Almost done!")
    
    run_command(options, "%s/run_tale_export.sh" % options["env"]["PBX_WORKFLOW_PATH"])
    
    print("Done!")
