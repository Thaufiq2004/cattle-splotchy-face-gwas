#!/usr/bin/env python3
"""Small glue steps between PLINK 2 commands (run from the repository root).

Sub-commands
------------
unique-ids   Give every variant a unique ID (chr:pos:a1:a2:row) -- 364 IDs are duplicated in the raw data.
phenotype    Convert the 0/1 phenotype file to PLINK coding (1 = control, 2 = case).
pcs          Turn a PLINK 2 .eigenvec file into a covariate file with columns PC1..PCn.
remove-list  Write the list of animals removed as a duplicated genotype (IID 15 and 23).
"""
import argparse
import shutil
import pandas as pd


def unique_ids(src_prefix, out_prefix):
    bim = pd.read_csv(f"{src_prefix}.bim", sep="\t", header=None, dtype=str)
    bim[1] = (bim[0] + ":" + bim[3] + ":" + bim[4] + ":" + bim[5] + ":"
              + (bim.index + 1).astype(str))
    assert bim[1].is_unique
    bim.to_csv(f"{out_prefix}.bim", sep="\t", header=False, index=False)
    shutil.copyfile(f"{src_prefix}.bed", f"{out_prefix}.bed")
    shutil.copyfile(f"{src_prefix}.fam", f"{out_prefix}.fam")
    print(f"wrote {out_prefix}.{{bed,bim,fam}} ({len(bim)} variants, all IDs unique)")


def phenotype(src, out):
    ph = pd.read_csv(src, sep=r"\s+")           # FID IID Splotchy (0 = control, 1 = case)
    assert set(ph["Splotchy"].unique()) <= {0, 1}
    ph["Splotchy"] = ph["Splotchy"] + 1         # PLINK: 1 = control, 2 = case
    ph.to_csv(out, sep="\t", index=False)
    print(f"wrote {out}: {int((ph.Splotchy == 2).sum())} cases, {int((ph.Splotchy == 1).sum())} controls")


def pcs(eigenvec, out, n):
    ev = pd.read_csv(eigenvec, sep=r"\s+", dtype=str)
    ev = ev.rename(columns={"#FID": "FID"})
    keep = ["FID", "IID"] + [f"PC{i}" for i in range(1, n + 1)]
    ev[keep].to_csv(out, sep="\t", index=False)
    print(f"wrote {out} ({len(ev)} animals, PC1-PC{n})")


def remove_list(out):
    with open(out, "w") as fh:
        fh.write("#FID\tIID\n0\t15\n0\t23\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("unique-ids"); a.add_argument("src_prefix"); a.add_argument("out_prefix")
    b = sub.add_parser("phenotype");  b.add_argument("src"); b.add_argument("out")
    c = sub.add_parser("pcs");        c.add_argument("eigenvec"); c.add_argument("out"); c.add_argument("--n", type=int, default=3)
    d = sub.add_parser("remove-list"); d.add_argument("out")
    args = ap.parse_args()
    if args.cmd == "unique-ids": unique_ids(args.src_prefix, args.out_prefix)
    elif args.cmd == "phenotype": phenotype(args.src, args.out)
    elif args.cmd == "pcs": pcs(args.eigenvec, args.out, args.n)
    else: remove_list(args.out)
