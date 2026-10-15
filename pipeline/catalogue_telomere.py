"""Telomere block caller for end windows in tip frame (contig tip = position 0, reading inward; a correctly oriented telomere
reads (CCCTAA)n). Same algorithm and parameters as the caller used to annotate the Release 2 tip table (junk_len, telo_len):
first position within `maxjunk` bases where 3 clean CCCTAA units start a 100 bp window >= 90% canonical; the block extends while
200 bp windows stay >= 60% canonical. junk_len = block start, telo_len = block end - start."""
import numpy as np
def _rots(m): return {m[i:] + m[:i] for i in range(len(m))}
CAN = _rots('CCCTAA')
def _covmask(s, motifs, k=6):
    m = np.zeros(len(s), bool)
    for i in range(len(s) - k + 1):
        if s[i:i + k] in motifs: m[i:i + k] = True
    return m
def _frac(cs, i, w, n):
    j = min(i + w, n); return (cs[j] - cs[i]) / max(1, j - i)
def find_block(s, motifs=CAN, W=100, clean_thr=0.9, ext_thr=0.6, maxjunk=8000):
    s = s.upper(); n = len(s); can = _covmask(s[:min(n, maxjunk + 60000)], motifs); n2 = len(can)
    cs = np.concatenate([[0], np.cumsum(can)]); start = None
    for i in range(0, min(maxjunk, n2 - W)):
        u = s[i:i + 6]
        if u in motifs and s[i:i + 18] == u * 3 and _frac(cs, i, W, n2) >= clean_thr: start = i; break
    if start is None: return None
    j = start
    while j + 200 <= n2 and _frac(cs, j, 200, n2) >= ext_thr: j += 50
    inner = min(j + 200, n2); idx = np.where(can[start:inner])[0]
    end = start + idx[-1] + 1 if len(idx) else start
    u = s[start:start + 6]; k = 0
    while k < 6 and start - k - 1 >= 0 and s[start - k - 1] == u[(-k - 1) % 6]: k += 1
    return dict(junk_len=start - k, telo_len=end - (start - k))
