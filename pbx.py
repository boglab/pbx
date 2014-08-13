import datetime
import subprocess
import os
import os.path
import ConfigParser
import shlex

def run_command(options, command):
    #subprocess.check_call(["bash", "-x", "-e"] + shlex.split(command))
    subprocess.check_call(["bash", "-e"] + shlex.split(command))

def validate_options(options):
    pass

def print_utils_and_citations():
    pass

def get_options(base_path):
    
    config = ConfigParser.RawConfigParser()
    config.read(os.path.join(base_path, "config.ini"))
    
    options = {}
    
    options["stop_after"] = "export"
    
    if config.has_option('General', 'stop_after'):
        options["stop_after"] = config.get('General', 'stop_after')
    
    options["dry_run"] = False
    
    if config.has_option('General', 'dry_run'):
        options["dry_run"] = config.getboolean('General', 'dry_run')
    
    # Environment
    
    options["env"] = {
        # Auto
        "PBX_BASE_PATH": base_path,
        "PBX_SCRIPTS_PATH": os.path.join(base_path, "scripts"),
        "PBX_PROTOCOL_TEMPLATES_PATH": os.path.join(base_path, "protocol_templates"),
        "PBX_WORKFLOW_PATH": os.path.join(base_path, "workflow"),
        # User config
        "PBX_SMRTANALYSIS_PATH": config.get('General', 'smrtanalysis_path'),
        "PBX_MUMMER_PATH": config.get('General', 'mummer_path'),
        "PBX_RAW_READS_PATH": config.get('General', 'raw_reads_path'),
        "PBX_TALE_SEQS_WHITELISTING": os.path.join(base_path, "tale_seqs", "whitelisting", config.get('Whitelisting', 'tale_seqs_file_whitelisting')),
        "PBX_TALE_SEQS_EXPORT": os.path.join(base_path, "tale_seqs", "exporter", config.get('Export', 'tale_seqs_file_export')),
        "PBX_TALE_SEQS_EXPORT_BOUNDARIES": os.path.join(base_path, "tale_seqs", "exporter", config.get('Export', 'tale_seqs_file_export_boundaries')),
    }
    
    if config.has_option('General', 'results_path'):
        options["env"]["PBX_RESULTS_PATH"] = config.get('General', 'results_path')
    else:
        options["env"]["PBX_RESULTS_PATH"] = os.path.join(base_path, "results", "results_" + datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S'))

    # Other auto set options
    
    options["env"]["PBX_RESULTS_PATH_WHITELISTING"] = os.path.join(options["env"]["PBX_RESULTS_PATH"], "blasr_raw_to_bls_tals_alignments")
    options["env"]["PBX_RESULTS_PATH_PREASSEMBLY"] = os.path.join(options["env"]["PBX_RESULTS_PATH"], "preassembly_job")
    options["env"]["PBX_RESULTS_PATH_ASSEMBLY"] = os.path.join(options["env"]["PBX_RESULTS_PATH"], "minimo_results")
    options["env"]["PBX_RESULTS_PATH_RESEQUENCING"] = os.path.join(options["env"]["PBX_RESULTS_PATH"], "resequencing")
    options["env"]["PBX_RESULTS_PATH_RESEQUENCED_CONTIG_ASSEMBLY"] = os.path.join(options["env"]["PBX_RESULTS_PATH"], "combine_resequenced_tals")
    
    return options

if __name__ == "__main__":
    
    base_path = os.path.dirname(os.path.realpath(__file__))
    
    options = get_options(base_path)
    
    os.environ.update(options["env"])
    
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
