import random, time
from typing import Set

'''from scripts.kmer_candidates_from_faa import encode_kmer_aa

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"  # 20 canonical AAs
AA_TO_INT = {ch: i for i, ch in enumerate(AA_ALPHABET)}
AA_TABLE = [255] * 256   # 255 = sentinel for "invalid AA"
for aa, i in AA_TO_INT.items():
    AA_TABLE[ord(aa)] = i

BASE = len(AA_ALPHABET)

def random_seq(n: int) -> str:
    return "".join(random.choice(AA_ALPHABET) for _ in range(n))

# ---- your original ----
def kmers_encoded_set_orig(seq: str, k: int) -> Set[int]:
    out: Set[int] = set()
    n = len(seq)
    for i in range(n - k + 1):
        code = encode_kmer_aa(seq[i:i+k])
        if code != -1:
            out.add(code)
    return out

# ---- optimized rolling+table (assuming AA_TABLE already built) ----
def kmers_encoded_set_fast(seq: str, k: int) -> Set[int]:
    n = len(seq)
    if k <= 0 or n < k:
        return set()

    s = seq.encode("ascii")
    out: Set[int] = set()
    out_add = out.add
    tbl = AA_TABLE
    base = BASE
    high = base ** (k - 1)

    code = 0
    valid_len = 0

    for i, b in enumerate(s):
        x = tbl[b]
        if x == 255:
            code = 0
            valid_len = 0
            continue

        if valid_len < k:
            code = code * base + x
            valid_len += 1
            if valid_len == k:
                out_add(code)
        else:
            old = tbl[s[i - k]]
            code = (code - old * high) * base + x
            out_add(code)

    return out

def bench(fn, seqs, k, reps=3):
    # warmup
    fn(seqs[0], k)
    t0 = time.perf_counter()
    for _ in range(reps):
        for s in seqs:
            fn(s, k)
    t1 = time.perf_counter()
    return (t1 - t0)

if __name__ == "__main__":
    random.seed(0)
    seqs = [random_seq(2000) for _ in range(200)]  # adjust sizes to match your data
    k = 5

    t_orig = bench(kmers_encoded_set_orig, seqs, k)
    t_fast = bench(kmers_encoded_set_fast, seqs, k)

    print(f"orig: {t_orig:.3f}s")
    print(f"fast: {t_fast:.3f}s")
    print(f"speedup: {t_orig / t_fast:.1f}x")'''
