"""Write results/MANIFEST.tsv: md5 and sha256 of every released file of this track (results/, figures/).
Usage: python scripts/checksums.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jcommon as J
rows = []
for d in (J.RES, J.FIG):
    for f in sorted(os.listdir(d)):
        p = os.path.join(d, f)
        if f == 'MANIFEST.tsv' or not os.path.isfile(p): continue
        rows.append((os.path.relpath(p, J.TRACK), os.path.getsize(p), J.md5(p), J.sha256(p)))
with open(os.path.join(J.RES, 'MANIFEST.tsv'), 'w') as f:
    f.write('file\tbytes\tmd5\tsha256\n')
    for r in rows: f.write('\t'.join(map(str, r)) + '\n')
print(len(rows), 'files')
