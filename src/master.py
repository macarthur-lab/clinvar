"""
Alternate implementation of master.bash with improved logging, skipping of commands if output files are up-to-date, parallelization, etc.

Run with -h to see all options.
"""

import configargparse
from datetime import datetime
import ftplib
import os
import sys

try:
    import pypez
    import pysam
    import pandas   # make sure all dependencies are installed
except ImportError as e:
    sys.exit("ERROR: Python module not installed. %s. Please run 'pip install -r requirements.txt' " % e)

p = configargparse.getArgParser()
g = p.add_argument_group('main args')
g.add("-R", "--reference-genome", help="b37 .fa genome reference file", required=True)
g.add("-E", "--exac-sites-vcf",  help="ExAC sites vcf file. If specified, a clinvar table with extra ExAC fields will also be created.")

pypez.init_command_line_args()
args = p.parse_args()

if not os.path.isfile(args.reference_genome):
    p.error("genome reference: file not found: %s" % args.reference_genome)
reference_genome = args.reference_genome

if args.exac_sites_vcf:
    if not os.path.isfile(args.exac_sites_vcf):
	p.error("ExAC sites vcf: file not found: %s" % args.exac_sites_vcf)
    if not os.path.isfile(args.exac_sites_vcf + ".tbi"):
        p.error("ExAC sites vcf: tabix index not found: %{s}.tbi" % args.exac_sites_vcf)

def get_remote_file_changed_time(ftp_host, ftp_path):
    """Returns time modified in seconds since the epoch"""

    try:
        ftp = ftplib.FTP(ftp_host)
        ftp.login()
        response = ftp.sendcmd("MDTM " + ftp_path)
        last_changed_time = datetime.strptime(response[4:], "%Y%m%d%H%M%S")
        return int(last_changed_time.strftime("%s"))  #.strftime("%d %B %Y %H:%M:%S")
    except Exception as e:
        print("ERROR: %s" % e)
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
print("Checking for new clinvar release")
download_if_changed(jr, "ClinVarFullRelease_00-latest.xml.gz",  "ftp.ncbi.nlm.nih.gov", "/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz")
download_if_changed(jr, "variant_summary.txt.gz",  "ftp.ncbi.nlm.nih.gov", "/pub/clinvar/tab_delimited/variant_summary.txt.gz")
jr.run()

job = pypez.Job()

# extract the GRCh37 coordinates, mutant allele, MeasureSet ID and PubMed IDs from it. This currently takes about 20 minutes.
job.add("python -u IN:parse_clinvar_xml.py -x IN:ClinVarFullRelease_00-latest.xml.gz -o OUT:clinvar_table_raw.tsv")

# normalize (convert to minimal representation and left-align)
# the normalization code is in a different repo (useful for more than just clinvar) so here I just wget it:
job.add("wget -N https://raw.githubusercontent.com/ericminikel/minimal_representation/master/normalize.py")
job.add("python -u normalize.py -R IN:%(reference_genome)s < IN:clinvar_table_raw.tsv > OUT:clinvar_table_normalized.tsv" % locals())

# join information from the tab-delimited summary to the normalized genomic coordinates
job.add("Rscript IN:join_data.R", input_filenames=['clinvar_table_normalized.tsv'], output_filenames=['clinvar_combined.tsv'])

# now sort again by genomic coordinates (because R's merge function ruins this)
job.add("(cat IN:clinvar_combined.tsv | head -1 > OUT:clinvar_combined_sorted.tsv ) && " + # header row
    "(cat IN:clinvar_combined.tsv | tail -n +2 | egrep -v \"^[XYM]\" | sort -k1,1n -k2,2n -k3,3 -k4,4 >> OUT:clinvar_combined_sorted.tsv ) && " + # numerically sort chroms 1-22
    "(cat IN:clinvar_combined.tsv | tail -n +2 | egrep \"^[XYM]\" | sort -k1,1 -k2,2n -k3,3 -k4,4 >> OUT:clinvar_combined_sorted.tsv )")     # lexicogaraphically sort non-numerical chroms at end

# now de-dup _again_, because the tab-delimited summary contains dups
job.add("python -u IN:dedup_clinvar.py < IN:clinvar_combined_sorted.tsv | tee clinvar.tsv | bgzip -c > OUT:clinvar.tsv.gz")  # clinvar_combined_sorted_dedup.tsv.gz
job.add("tabix -S 1 -s 1 -b 2 -e 2 IN:clinvar.tsv.gz", output_filenames=["clinvar.tsv.gz.tbi"])
job.add("cp IN:clinvar.tsv.gz IN:clinvar.tsv.gz.tbi ../output", output_filenames=["../output/clinvar.tsv", "../output/clinvar.tsv.gz", "../output/clinvar.tsv.gz.tbi"])

# create vcf
job.add("python -u IN:clinvar_table_to_vcf.py IN:clinvar.tsv | bgzip -c > OUT:clinvar.vcf.gz")  # create compressed version
job.add("tabix IN:clinvar.vcf.gz", output_filenames=["clinvar.vcf.gz.tbi"])
job.add("cp IN:clinvar.vcf.gz IN:clinvar.vcf.gz.tbi ../output", output_filenames=["../output/clinvar.vcf.gz", "../output/clinvar.vcf.gz.tbi"])

# create tsv table with extra fields from ExAC: filter, ac_adj, an_adj, popmax_ac, popmax_an, popmax
if args.exac_sites_vcf:
    normalized_vcf = os.path.basename(args.exac_sites_vcf).split('.vcf')[0] + ".normalized.vcf.gz"
    job.add("vt decompose -s IN:%s | vt normalize -r IN:%s - | bgzip -c > OUT:%s" % (args.exac_sites_vcf, args.reference_genome, normalized_vcf))
    job.add("tabix IN:"+normalized_vcf, output_filenames=[normalized_vcf+".tbi"])
    job.add("python -u IN:add_exac_fields.py -i IN:clinvar.tsv -e IN:%(normalized_vcf)s | bgzip -c > OUT:clinvar_with_exac.tsv.gz" % locals())
    job.add("tabix -S 1 -s 1 -b 2 -e 2 IN:clinvar_with_exac.tsv.gz", output_filenames=["clinvar_with_exac.tsv.gz.tbi"])
    job.add("cp IN:clinvar_with_exac.tsv.gz IN:clinvar_with_exac.tsv.gz.tbi ../output", output_filenames=["../output/clinvar_with_exac.tsv.gz", "../output/clinvar_with_exac.tsv.gz.tbi"])

# create uncompressed example files that contain the 1st 5000 lines of the compressed tsvs so people can easily see typical values online on github
job.add("gunzip -c IN:clinvar.vcf.gz | head -n 5000 > OUT:../output/clinvar_example_5000_rows.vcf")
job.add("gunzip -c IN:clinvar.tsv.gz | head -n 5000 > OUT:../output/clinvar_example_5000_rows.tsv")
job.add("gunzip -c IN:clinvar_with_exac.tsv.gz | head -n 5000 > OUT:../output/clinvar_with_exac_example_5000_rows.tsv")

# create a stats file that summarizes some of the columns
# chrom      pos        ref        alt        mut        measureset_id      symbol     clinical_significance      review_status      hgvs_c     hgvs_p     all_submitters     all_traits         all_pmids          inheritance_modes          age_of_onset        prevalence         disease_mechanism          origin     xrefs      gold_stars         pathogenic         conflicted
job.add("""echo Total: $(gunzip -c clinvar.tsv.gz | tail -n +2 | wc -l) > OUT:clinvar_stats.txt &&
for i in 8 9 16 17 18 19 21; do 
    echo ================ >> OUT:clinvar_stats.txt ;
    gunzip -c IN:clinvar.tsv.gz | head -n 1 | cut -f $i >> OUT:clinvar_stats.txt ;
    gunzip -c IN:clinvar.tsv.gz | tail -n +2 | cut -f $i | tr ';' '\n' | sort | uniq -c | sort -r -n >> OUT:clinvar_stats.txt ;
done
""", input_filenames=["clinvar.tsv.gz"])

job.add("cp IN:clinvar_stats.txt OUT:../output/clinvar_stats.txt")

# run the above commands
jr.run(job)
