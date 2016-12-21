'''
Created on 21 Dec 2016

@author: xzhang13 AT ic.ac.uk
'''
"""
Download the latest clinvar xml and parse it to table, then normalize and sort. The table records variant-condition pair info. 

Run with -h to see all options.
"""

import configargparse
from datetime import datetime
import ftplib
import os
import sys
from distutils import spawn

try:
    import pypez
    import pysam
    import pandas   # make sure all dependencies are installed
except ImportError as e:
    sys.exit("ERROR: Python module not installed. %s. Please run 'pip install -r requirements.txt' " % e)
for executable in ['wget', 'Rscript', 'tabix', 'vt']:
    assert spawn.find_executable(executable), "Command %s not found, see README" % executable

p = configargparse.getArgParser()
g = p.add_argument_group('main args')
g.add("-R", "--reference-genome", help="b37 .fa genome reference file", required=True)
#g.add("-E", "--exac-sites-vcf",  help="ExAC sites vcf file. If specified, a clinvar table with extra ExAC fields will also be created.")
g.add("-X", "--clinvar-xml", help="The local filename of the ClinVarFullRelase.xml.gz file. If not set, grab the latest from NCBI.")
#g.add("-S", "--clinvar-variant-summary-table", help="The local filename of the variant_summary.txt.gz file. If not set, grab the latest from NCBI.")

pypez.init_command_line_args()
args = p.parse_args()

if not os.path.isfile(args.reference_genome):
    p.error("genome reference: file not found: %s" % args.reference_genome)
reference_genome = args.reference_genome


def get_remote_file_changed_time(ftp_host, ftp_path):
    """Returns time modified in seconds since the epoch"""
    print("Retrieving last-changed time for %s" % os.path.join(ftp_host, ftp_path))
    try:
        ftp = ftplib.FTP(ftp_host, timeout=3)  # timeout = 3 seconds
        ftp.login()
        response = ftp.sendcmd("MDTM " + ftp_path)
        last_changed_time = datetime.strptime(response[4:], "%Y%m%d%H%M%S")
        return int(last_changed_time.strftime("%s"))  #.strftime("%d %B %Y %H:%M:%S")
    except Exception as e:
        print("ERROR: retrieving last-changed time for %s: %s" % (os.path.join(ftp_host, ftp_path), e))
        return 0


def download_if_changed(job_runner, local_path, ftp_host, ftp_path):
    remote_changed_time = get_remote_file_changed_time(ftp_host, ftp_path)
    local_changed_time = os.path.getmtime(local_path) if os.path.isfile(local_path) else 0
    ftp_address = "ftp://%s/%s" % (ftp_host, ftp_path)
    if remote_changed_time > local_changed_time:
        job_runner.add_parallel(pypez.Job("wget %s -O OUT:%s" % (ftp_address, local_path)))
    else:
        print("Local copy of %s is up to date. The remote version hasn't changed since %s" % (ftp_address, datetime.fromtimestamp(remote_changed_time)))
    #ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz

jr = pypez.JobRunner()

if args.clinvar_xml:
    if not os.path.isfile(args.clinvar_xml):
        p.error("ClinVar XML specified but not found: %s" % args.clinvar_xml)
    if not args.clinvar_xml.endswith('.gz'):
        p.error("ClinVar XML expected to be gzipped: %s" % args.clinvar_xml)
    clinvar_xml = args.clinvar_xml
else:
    print("Checking for new clinvar release")
    clinvar_xml = "ClinVarFullRelease_00-latest.xml.gz"
    download_if_changed(jr, clinvar_xml,  "ftp.ncbi.nlm.nih.gov", "/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz")

jr.run()

job = pypez.Job()

# extract the GRCh37 coordinates, mutant allele, MeasureSet ID and PubMed IDs from it. This currently takes about 20 minutes.
job.add("python -u IN:parse_clinvar_xml.py -x IN:%s -o OUT:clinvar_table_raw.tsv -e parse_clinvar_table_raw.logs" % clinvar_xml)

# normalize (convert to minimal representation and left-align)
# the normalization code is in a different repo (useful for more than just clinvar) so here I just wget it:
job.add("wget -N https://raw.githubusercontent.com/ericminikel/minimal_representation/master/normalize.py")
job.add("python -u normalize.py -R IN:%(reference_genome)s < IN:clinvar_table_raw.tsv > OUT:clinvar_table_normalized.tsv" % locals())


# sort 
job.add("ex -s +'bufdo!v/\S/d' -cxa clinvar_table_normalized.tsv") #remove the empty lines
job.add("(cat IN:clinvar_table_normalized.tsv | head -1 > OUT:clinvar_sorted.tsv ) && " + # header row
    "(cat IN:clinvar_table_normalized.tsv | tail -n +2 | egrep -v \"^[XYM]\" | sort -k1,1n -k2,2n -k3,3 -k4,4 >> OUT:clinvar_sorted.tsv ) && " + # numerically sort chroms 1-22
    "(cat IN:clinvar_table_normalized.tsv | tail -n +2 | egrep \"^[XYM]\" | sort -k1,1 -k2,2n -k3,3 -k4,4 >> OUT:clinvar_sorted.tsv )")     # lexicogaraphically sort non-numerical chroms at end
job.add("cp clinvar_sorted.tsv ../clinvar.tsv")



# run the above commands
jr.run(job)
