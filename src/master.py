"""
Alternate implementation of master.bash with improved logging, skipping of commands if output files are up-to-date, parallelization, etc.

Run with -h to see all options.
"""

import configargparse
from datetime import datetime
import ftplib
import os
import sys
from distutils import spawn

try:
    import configargparse
    import pypez
    import pysam
    import pandas   # make sure all dependencies are installed
except ImportError as e:
    sys.exit("ERROR: Python module not installed. %s. Please run 'pip install -r requirements.txt' " % e)
for executable in ['wget', 'tabix', 'vt']:
    assert spawn.find_executable(executable), "Command %s not found, see README" % executable

p = configargparse.getArgParser()
g = p.add_argument_group('main args')
g.add("--b37-genome", help="b37 .fa genome reference file", default=None, required=False)
g.add("--b38-genome", help="b38 .fa genome reference file. NOTE: chromosome names must be like '1', '2'.. 'X', 'Y', 'MT'.", default=None, required=False)
g.add("-X", "--clinvar-xml", help="The local filename of the ClinVarFullRelase.xml.gz file. If not set, grab the latest from NCBI.")
g.add("-S", "--clinvar-variant-summary-table", help="The local filename of the variant_summary.txt.gz file. If not set, grab the latest from NCBI.")
g.add("-E", "--exac-sites-vcf",  help="ExAC sites vcf file. If specified, a clinvar table with extra ExAC fields will also be created.")
g.add("-GE", "--gnomad-exome-sites-vcf",  help="gnomAD exome sites vcf file. If specified, a clinvar table with extra gnomAD exome info fields will also be created.")
g.add("-GG", "--gnomad-genome-sites-vcf",  help="gnomAD genome sites vcf file. If specified, a clinvar table with extra gnomAD genome info fields will also be created.")
g.add("--output-prefix", default="../output/", help="Final output files will have this prefix")
g.add("--tmp-dir", default="./output_tmp", help="Temporary output files will have this prefix")
g = p.add_mutually_exclusive_group()
g.add("--single-only", dest="single_or_multi", action="store_const", const="single", help="Only generate the single-variant tables")
g.add("--multi-only", dest="single_or_multi", action="store_const", const="multi", help="Only generate the multi-variant tables")


pypez.init_command_line_args()
args = p.parse_args()
for key, value in args.__dict__.items():
    print("%s=%s" % (key, value))

reference_genomes = {'b37': args.b37_genome, 'b38': args.b38_genome}
clinvar_xml = args.clinvar_xml
#if clinvar_xml and not os.path.isfile(clinvar_xml)
exac_sites_vcf = args.exac_sites_vcf
gnomad_exome_sites_vcf = args.gnomad_exome_sites_vcf
gnomad_genome_sites_vcf = args.gnomad_genome_sites_vcf
clinvar_variant_summary_table = args.clinvar_variant_summary_table
output_prefix = args.output_prefix

tmp_dir = args.tmp_dir
os.system("mkdir -p " + tmp_dir)

if reference_genomes['b37'] is None and reference_genomes['b38'] is None:
    p.error("At least one genome reference file is required")

for key, path in reference_genomes.items():
    if path is not None and not os.path.isfile(path):
        p.error("%s genome reference: file not found: %s" % (key, path))

for label, vcf_path in (('gnomad_genomes', gnomad_genome_sites_vcf), ('gnomad_exomes', gnomad_exome_sites_vcf), ('exac_v1', exac_sites_vcf)):
    if not vcf_path:
        continue

    if not os.path.isfile(vcf_path):
        p.error(label+" sites vcf: file not found: %s" % vcf_path)
    if not os.path.isfile(vcf_path + ".tbi"):
        p.error(label+" sites vcf: tabix index not found: %s.tbi" % vcf_path)


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
        print("Local copy of %s is out of date. The remote version changed on %s" % (ftp_address, datetime.fromtimestamp(remote_changed_time)))
        job_runner.add_parallel(pypez.Job("wget %s -O OUT:%s" % (ftp_address, local_path)))
    else:
        print("Local copy of %s is up to date. The remote version hasn't changed since %s" % (ftp_address, datetime.fromtimestamp(remote_changed_time)))
        #ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz

jr = pypez.JobRunner()

if clinvar_xml:
    if not os.path.isfile(clinvar_xml):
        p.error("ClinVar XML specified but not found: %s" % clinvar_xml)
    if not clinvar_xml.endswith('.gz'):
        p.error("ClinVar XML expected to be gzipped: %s" % clinvar_xml)
    clinvar_xml = clinvar_xml
else:
    print("Checking for new clinvar release")
    clinvar_xml = "ClinVarFullRelease_00-latest.xml.gz"
    download_if_changed(jr, clinvar_xml,  "ftp.ncbi.nlm.nih.gov", "/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz")

if clinvar_variant_summary_table:
    if not os.path.isfile(clinvar_variant_summary_table):
        p.error("ClinVar variant summary table specified but not found: %s" % clinvar_variant_summary_table)
    if not clinvar_variant_summary_table.endswith('.gz'):
        p.error("ClinVar variant summary table expected to be gzipped: %s" % clinvar_variant_summary_table)
    variant_summary_table = clinvar_variant_summary_table
else:
    print("Checking for new clinvar release")
    variant_summary_table = "variant_summary.txt.gz"
    download_if_changed(jr, variant_summary_table,  "ftp.ncbi.nlm.nih.gov", "/pub/clinvar/tab_delimited/variant_summary.txt.gz")

jr.run()

job = pypez.Job()

# normalize (convert to minimal representation and left-align)
# the normalization code is in a different repo (useful for more than just clinvar) so here I just wget it:
job.add("wget -N https://raw.githubusercontent.com/ericminikel/minimal_representation/master/normalize.py")

for genome_build in ('b37', 'b38'):
    # extract the GRCh37 coordinates, mutant allele, MeasureSet ID and PubMed IDs from it. This currently takes about 20 minutes.
    genome_build_id = genome_build.replace('b', 'GRCh')
    reference_genome = reference_genomes[genome_build]
    if reference_genome is None:
        print("Skippping steps to generate %s tables since reference genome not given." % genome_build)
        continue

    job.add(("python -u IN:parse_clinvar_xml.py "
            "-x IN:%(clinvar_xml)s "
            "-g %(genome_build_id)s "
            "-o OUT:%(tmp_dir)s/clinvar_table_raw.single.%(genome_build)s.tsv "
            "-m OUT:%(tmp_dir)s/clinvar_table_raw.multi.%(genome_build)s.tsv") % locals())

    for is_multi in (True, False):  # multi = clinvar submission that describes multiple alleles (eg. compound het, haplotypes, etc.)
        single_or_multi = 'multi' if is_multi else 'single'
        if args.single_or_multi and single_or_multi != args.single_or_multi:
            print("Skippping steps to generate %s table." % single_or_multi)
            continue
            
        fsuffix = "%(single_or_multi)s.%(genome_build)s" % locals() # file suffix
        output_dir = '%(output_prefix)s%(genome_build)s/%(single_or_multi)s' % locals()
        os.system('mkdir -p ' + output_dir)

        # normalize variants  (use grep -v '^$' to remove empty rows)
        job.add("python -u normalize.py -R IN:%(reference_genome)s < IN:%(tmp_dir)s/clinvar_table_raw.%(fsuffix)s.tsv | grep -v ^$ | bgzip -c > OUT:%(tmp_dir)s/clinvar_table_normalized.%(fsuffix)s.tsv.gz" % locals())

        # sort
        job.add(("cat " +
            "<(gunzip -c IN:%(tmp_dir)s/clinvar_table_normalized.%(fsuffix)s.tsv.gz | head -1) "  # header row
            "<(gunzip -c IN:%(tmp_dir)s/clinvar_table_normalized.%(fsuffix)s.tsv.gz | tail -n +2 | egrep -v \"^[XYM]\" | sort -k1,1n -k2,2n -k3,3 -k4,4 ) " + # numerically sort chroms 1-22
            "<(gunzip -c IN:%(tmp_dir)s/clinvar_table_normalized.%(fsuffix)s.tsv.gz | tail -n +2 | egrep \"^[XYM]\" | sort -k1,1 -k2,2n -k3,3 -k4,4 ) " +  #sort chroms X,Y,M 
            " | bgzip -c > OUT:%(tmp_dir)s/clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz") % locals())   # lexicogaraphically sort non-numerical chroms at end

        # tabix and copy to output dir
        job.add("tabix -S 1 -s 1 -b 2 -e 2 IN:%(tmp_dir)s/clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz" % locals(), output_filenames=["%(tmp_dir)s/clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz.tbi" % locals()])
        job.add("cp IN:%(tmp_dir)s/clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz IN:%(tmp_dir)s/clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz.tbi %(output_dir)s" % locals(), output_filenames=[
            "%(output_dir)s/clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz" % locals(),
            "%(output_dir)s/clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz.tbi" % locals()
            ])

        # group by allele, since clinvar_allele_trait_pairs.*.tsv will have more than 1 record for some alleles
        job.add("python -u IN:group_by_allele.py -i IN:%(tmp_dir)s/clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz | bgzip -c > OUT:%(tmp_dir)s/clinvar_alleles_grouped.%(fsuffix)s.tsv.gz" % locals())

        # join information from the tab-delimited summary to the normalized genomic coordinates
        job.add("python IN:join_variant_summary_with_clinvar_alleles.py "
                "IN:%(variant_summary_table)s "
                "IN:%(tmp_dir)s/clinvar_alleles_grouped.%(fsuffix)s.tsv.gz "
                "OUT:%(tmp_dir)s/clinvar_alleles_combined.%(fsuffix)s.tsv.gz "
                "%(genome_build_id)s" % locals())

        # sort again by genomic coordinates
        job.add(("cat " +
            "<(gunzip -c IN:%(tmp_dir)s/clinvar_alleles_combined.%(fsuffix)s.tsv.gz | head -1) "  # header row
            "<(gunzip -c IN:%(tmp_dir)s/clinvar_alleles_combined.%(fsuffix)s.tsv.gz | tail -n +2 | egrep -v \"^[XYM]\" | sort -k1,1n -k2,2n -k3,3 -k4,4 ) " + # numerically sort chroms 1-22
            "<(gunzip -c IN:%(tmp_dir)s/clinvar_alleles_combined.%(fsuffix)s.tsv.gz | tail -n +2 | egrep \"^[XYM]\" | sort -k1,1 -k2,2n -k3,3 -k4,4 ) " +  #sort chroms X,Y,M 
            " | bgzip -c > OUT:%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.tsv.gz") % locals())   # lexicogaraphically sort non-numerical chroms at end

        # tabix and copy to output dir
        job.add("tabix -S 1 -s 1 -b 2 -e 2 IN:%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.tsv.gz" % locals(), output_filenames=["%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.tsv.gz.tbi" % locals()])
        job.add("cp IN:%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.tsv.gz IN:%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.tsv.gz.tbi %(output_dir)s/"  % locals(),
                output_filenames=[
                    "%(output_dir)s/clinvar_alleles.%(fsuffix)s.tsv.gz" % locals(),
                    "%(output_dir)s/clinvar_alleles.%(fsuffix)s.tsv.gz.tbi" % locals()
                ])

        # create vcf
        job.add(("python -u IN:clinvar_table_to_vcf.py IN:%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.tsv.gz IN:%(reference_genome)s | bgzip -c > OUT:%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.vcf.gz") % locals())  # create compressed version
        job.add("tabix IN:%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.vcf.gz" % locals(), output_filenames=["%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.vcf.gz.tbi" % locals()])
        job.add("cp IN:%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.vcf.gz IN:%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.vcf.gz.tbi %(output_dir)s/" % locals(), output_filenames=[
            "%(output_dir)s/clinvar_alleles.%(fsuffix)s.vcf.gz" % locals(),
            "%(output_dir)s/clinvar_alleles.%(fsuffix)s.vcf.gz.tbi" % locals()])

        # create uncompressed example files that contain the 1st 750 lines of the compressed tsvs so people can easily see typical values online on github
        job.add("gunzip -c IN:%(tmp_dir)s/clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz | head -n 750 > OUT:%(output_dir)s/clinvar_allele_trait_pairs_example_750_rows.%(fsuffix)s.tsv" % locals())
        job.add("gunzip -c IN:%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.tsv.gz | head -n 750 > OUT:%(output_dir)s/clinvar_alleles_example_750_rows.%(fsuffix)s.tsv" % locals())
        job.add("gunzip -c IN:%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.vcf.gz | head -n 750 > OUT:%(output_dir)s/clinvar_alleles_example_750_rows.%(fsuffix)s.vcf" % locals())

        # create tsv table with extra fields from ExAC: filter, ac_adj, an_adj, popmax_ac, popmax_an, popmax
        if genome_build == "b37":
            for label, vcf_arg, vcf_path in (('gnomad_genomes', '-gg', gnomad_genome_sites_vcf), ('gnomad_exomes', '-ge', gnomad_exome_sites_vcf), ('exac_v1', '-e', exac_sites_vcf)):
                if not vcf_path:
                    continue
                script_name = "add_exac_fields.py" if label == "exac_v1" else "add_gnomad_fields.py"
                normalized_vcf = os.path.basename(vcf_path).split('.vcf')[0] + ".normalized.vcf.gz" % locals()
                job.add(("vt decompose -s IN:%(vcf_path)s | "
                         "vt normalize -r IN:%(reference_genome)s - | "
                         "bgzip -c > OUT:%(tmp_dir)s/%(normalized_vcf)s") % locals())
                job.add("tabix IN:%(tmp_dir)s/%(normalized_vcf)s" % locals(), output_filenames=["%(tmp_dir)s/%(normalized_vcf)s.tbi" % locals()])
                job.add(("python -u IN:%(script_name)s -i IN:%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.tsv.gz %(vcf_arg)s IN:%(tmp_dir)s/%(normalized_vcf)s | "
                         "bgzip -c > OUT:%(tmp_dir)s/clinvar_alleles_with_%(label)s.%(fsuffix)s.tsv.gz") % locals())
                job.add("tabix -S 1 -s 1 -b 2 -e 2 IN:%(tmp_dir)s/clinvar_alleles_with_%(label)s.%(fsuffix)s.tsv.gz" % locals(), output_filenames=["%(tmp_dir)s/clinvar_alleles_with_%(label)s.%(fsuffix)s.tsv.gz.tbi" % locals()])
                job.add("cp IN:%(tmp_dir)s/clinvar_alleles_with_%(label)s.%(fsuffix)s.tsv.gz IN:%(tmp_dir)s/clinvar_alleles_with_%(label)s.%(fsuffix)s.tsv.gz.tbi %(output_dir)s/" % locals(), output_filenames=[
                    "%(output_dir)s/clinvar_alleles_with_%(label)s.%(fsuffix)s.tsv.gz" % locals(),
                    "%(output_dir)s/clinvar_alleles_with_%(label)s.%(fsuffix)s.tsv.gz.tbi" % locals()])

                job.add("gunzip -c IN:%(tmp_dir)s/clinvar_alleles_with_%(label)s.%(fsuffix)s.tsv.gz | head -n 750 > OUT:%(output_dir)s/clinvar_alleles_with_%(label)s_example_750_rows.%(fsuffix)s.tsv" % locals())

        job.add(
            "python clinvar_alleles_stats.py "
            "IN:%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.tsv.gz "
            "> OUT:%(tmp_dir)s/clinvar_alleles_stats.%(fsuffix)s.txt" %
            locals(),
            input_filenames=[
                "%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.tsv.gz" % locals(),
                "clinvar_alleles_stats.py"])

        job.add("cp IN:%(tmp_dir)s/clinvar_alleles_stats.%(fsuffix)s.txt OUT:%(output_dir)s/clinvar_alleles_stats.%(fsuffix)s.txt" % locals())

        # run basic checks
        job.add("python IN:check_allele_table.py IN:%(tmp_dir)s/clinvar_alleles.%(fsuffix)s.tsv.gz" % locals())

# run the above commands
jr.run(job)
