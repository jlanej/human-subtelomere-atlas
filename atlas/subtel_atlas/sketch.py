"""Deterministic FracMinHash sketches of DNA (canonical k-mers, k <= 31, 2-bit packed, splitmix64 hash; kept if hash < 2^64/scale)
and Mash-style ANI between sketches: ANI = 1 + ln(2J / (1 + J)) / k, J = Jaccard of the two sketches."""
import numpy as np
_LUT = np.full(256, 4, np.uint8)
for i, c in enumerate(b'ACGT'): _LUT[c] = i
for i, c in enumerate(b'acgt'): _LUT[c] = i
def _splitmix(x):
    x = (x + np.uint64(0x9E3779B97F4A7C15))
    x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return x ^ (x >> np.uint64(31))
def sketch(seq, k=21, scale=10):
    a = _LUT[np.frombuffer(seq.encode(), np.uint8)]; n = len(a) - k + 1
    if n <= 0: return np.zeros(0, np.uint64)
    bad = (a == 4).astype(np.int32); badk = np.convolve(bad, np.ones(k, np.int32), 'valid') > 0
    v = np.where(a == 4, 0, a).astype(np.uint64); r = (np.uint64(3) - v)
    fw = np.zeros(n, np.uint64); rv = np.zeros(n, np.uint64)
    with np.errstate(over='ignore'):
        for j in range(k):
            fw = (fw << np.uint64(2)) | v[j:j + n]
            rv = rv | (r[j:j + n] << np.uint64(2 * j))
        h = _splitmix(np.minimum(fw, rv))
    keep = (~badk) & ((h < np.uint64(((1 << 64) - 1) // scale)) if scale > 1 else True)
    return np.unique(h[keep])
def ani_from_jaccard(J, k=21):
    J = np.asarray(J, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        a = 1 + np.log(2 * J / (1 + J)) / k
    return np.where(J > 0, a, 0.0)
