#!/usr/bin/env python
"""Rebuilds the perturb-seq knockdown signature matrices used by the STAT4 vignette panel d.

Opt-in: not part of make tables; run via make knn-signatures. Writes data/perturbseq/knn/signatures_<COND>.npy (float32, n_perturbed x 10282), signatures_<COND>.genes.txt, and var_gene_names.txt.

See DATA_AVAILABILITY.md, Perturb-seq Knockdown Signatures, for source, filtering, and provenance.
"""
import os, sys, time, urllib.request
import numpy as np, h5py

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
assert os.path.exists(os.path.join(ROOT, "Makefile")), f"repo root not found at {ROOT}"

OUT = os.path.join(ROOT, "data", "perturbseq", "knn")
S3  = ("https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/"
       "marson2025_data/GWCD4i.DE_stats.h5ad")
CONDS = ["Rest", "Stim8hr", "Stim48hr"]
MIN_CELLS = 50


class HTTPRangeFile:
    """Minimal seekable read-only file over HTTP range requests (for h5py header reads)."""
    def __init__(self, url, block=4*1024*1024):
        self.url=url; self.block=block; self.pos=0; self.cache={}
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60) as r:
            self.size=int(r.headers["Content-Length"])
    def _fetch(self, b):
        if b in self.cache: return self.cache[b]
        s=b*self.block; e=min(s+self.block, self.size)-1
        with urllib.request.urlopen(urllib.request.Request(self.url, headers={"Range":f"bytes={s}-{e}"}), timeout=120) as r:
            data=r.read()
        if len(self.cache)>256: self.cache.clear()
        self.cache[b]=data; return data
    def seek(self, off, whence=0):
        self.pos=(off if whence==0 else self.pos+off if whence==1 else self.size+off); return self.pos
    def tell(self): return self.pos
    def seekable(self): return True
    def readable(self): return True
    def read(self, n=-1):
        if n<0: n=self.size-self.pos
        end=min(self.pos+n, self.size); out=bytearray()
        for b in range(self.pos//self.block, (end-1)//self.block+1):
            data=self._fetch(b); s=b*self.block
            out+=data[max(self.pos,s)-s: min(end,s+len(data))-s]
        self.pos=end; return bytes(out)
    def readinto(self, b):
        d=self.read(len(b)); b[:len(d)]=d; return len(d)


def read_cat(hf, grp):
    o=hf[grp]
    if isinstance(o, h5py.Group) and "categories" in o:
        return np.asarray(o["categories"].asstr()[:])[o["codes"][:]]
    return o.asstr()[:] if o.dtype.kind in "OS" else o[:]


def main():
    os.makedirs(OUT, exist_ok=True)
    hf=h5py.File(HTTPRangeFile(S3), "r")
    lay=hf["layers"]["zscore"]; nrow, ncol=lay.shape; off=lay.id.get_offset()
    assert lay.chunks is None and off is not None, "zscore layer must be contiguous/uncompressed"
    var_name=np.asarray(hf["var/gene_name"].asstr()[:])
    cond=read_cat(hf, "obs/culture_condition")
    pert=read_cat(hf, "obs/target_contrast_gene_name")
    ont =read_cat(hf, "obs/ontarget_significant")
    ncell=np.asarray(hf["obs/n_cells_target"][:], float)
    ont_b=ont if ont.dtype==bool else np.array([str(x).lower() in ("true","1","1.0") for x in ont])

    np.savetxt(os.path.join(OUT, "var_gene_names.txt"), var_name, fmt="%s")

    # Staged via a float32 memmap since the source layer is float64 (~2.8GB); temp file removed only on success.
    tmp=np.memmap(os.path.join(OUT, "_zscore_tmp.dat"), dtype=np.float32, mode="w+", shape=(nrow, ncol))
    step=400; t0=time.time()
    for r0 in range(0, nrow, step):
        r1=min(r0+step, nrow); b0, b1=off+r0*ncol*8, off+r1*ncol*8-1
        for att in range(5):
            try:
                with urllib.request.urlopen(urllib.request.Request(S3, headers={"Range":f"bytes={b0}-{b1}"}), timeout=180) as rr:
                    buf=rr.read()
                break
            except Exception:
                # Retries S3 stream failures up to 5 times with linear backoff, then raises rather than write a truncated block.
                if att==4: raise
                time.sleep(2*(att+1))
        tmp[r0:r1,:]=np.frombuffer(buf, dtype="<f8").astype(np.float32).reshape(r1-r0, ncol)
        if r0 % 8000 < step: print(f"[stream] {r0}/{nrow} {time.time()-t0:.0f}s", flush=True)
    tmp.flush()

    for C in CONDS:
        sel=(cond==C) & ont_b & (ncell >= MIN_CELLS)
        ridx=np.where(sel)[0]; rg=pert[ridx]
        # return_index gives first-occurrence positions; re-sorting them preserves source row order, matching .genes.txt.
        _, keep=np.unique(rg, return_index=True); ridx=ridx[np.sort(keep)]; rg=pert[ridx]
        X=np.asarray(tmp[ridx,:], dtype=np.float32)
        np.save(os.path.join(OUT, f"signatures_{C}.npy"), X)
        np.savetxt(os.path.join(OUT, f"signatures_{C}.genes.txt"), rg, fmt="%s")
        print(f"[{C}] {X.shape} -> signatures_{C}.npy", flush=True)

    del tmp
    os.remove(os.path.join(OUT, "_zscore_tmp.dat"))
    print("DONE -- now re-run src/analysis/vignette/build_vignette_panelD_knn.py to re-derive panel d")


if __name__ == "__main__":
    main()
