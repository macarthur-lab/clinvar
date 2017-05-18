import sys
import re
import gzip

"""
Helper script to grab some measurets by their ID from the master XML for
testing purposes.
Usage:
    python grab_interesting_measuresets.py \
        <ClinVarFullRelease.xml.gz> \
        <comma-separated list of measureset IDs> \
        <out.xml.gz>
"""

in_xml = sys.argv[1]  # e.g. ClinVarFullRelease_2017-02.xml.gz
interesting_measuresets = set(sys.argv[2].split(","))
# ^ comma-separated list of interesting measuret IDs, e.g. 187175,188901
out_xml = sys.argv[3]  # where to write, e.g. interesting.xml.gz

measureset_id_regex = re.compile(r'ID="(\d+)"')

# input file could be gzipped or not, output file will have same status
if in_xml.endswith(".gz"):
    in_f = gzip.open(in_xml)
    if not out_xml.endswith(".gz"):
        out_xml += ".gz"
    out_f = gzip.open(out_xml, 'w')
else:
    in_f = open(in_xml)
    assert not out_xml.endswith(".gz")
    out_f = open(out_xml, 'w')

in_clinvarset = False
interesting = False
clinvarset = []
out_f.write(next(in_f))  # <?xml>
out_f.write(next(in_f))  # <RelaseSet>
out_f.write("\n")

for line in in_f:
    if line.startswith("<ClinVarSet"):
        in_clinvarset = True
    elif line.startswith("</ClinVarSet>"):
        if interesting:
            out_f.write("".join(clinvarset))
            out_f.write(line)
            out_f.write("\n")
        clinvarset = []
        in_clinvarset = False
        interesting = False
        continue
    else:
        if line.startswith("    <MeasureSet"):
            m = measureset_id_regex.search(line)
            interesting = (interesting or
                           (m and m.group(1) in interesting_measuresets))
    if in_clinvarset:
        clinvarset.append(line)

out_f.write("</ReleaseSet>\n")

in_f.close()
out_f.close()

