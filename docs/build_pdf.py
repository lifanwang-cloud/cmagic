#!/usr/bin/env python
"""Build docs/CMAGIC_COOKBOOK.pdf from the cookbook markdown (pandoc + xelatex).
Preprocessing maps unicode sub/superscripts and math symbols that common text
fonts lack into LaTeX math, leaving fenced code blocks verbatim (script-E -> Ecal).
Usage: python docs/build_pdf.py [input.md]  (default: the cookbook markdown)
"""
import re, subprocess, sys, os, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
src = sys.argv[1] if len(sys.argv) > 1 else (
    os.path.join(HERE, "CMAGIC_COOKBOOK.md") if os.path.exists(os.path.join(HERE, "CMAGIC_COOKBOOK.md"))
    else os.path.join(HERE, "METHOD.md"))
t = open(src).read()
sub = dict(zip("₀₁₂₃₄₅₆₇₈₉", "0123456789"))
parts = re.split(r"(```[\s\S]*?```)", t)
for i, p in enumerate(parts):
    if p.startswith("```"):
        parts[i] = p.replace("ℰ", "Ecal")
        continue
    # \ensuremath, not $..$ or \(..\): a $-span abutting a digit fails pandoc's
    # math parsing (the command then runs in TEXT mode -> missing text glyph),
    # and \( is swallowed by pandoc's backslash-escape rule. \ensuremath{..}
    # passes through as raw LaTeX and is safe in any mode and any adjacency.
    p = re.sub("[₀-₉]+", lambda m: r"\ensuremath{_{" + "".join(sub[c] for c in m.group(0)) + "}}", p)
    for a, b in [("≈", r"\\ensuremath{\\approx}"), ("≤", r"\\ensuremath{\\le}"),
                 ("≥", r"\\ensuremath{\\ge}"), ("≳", r"\\ensuremath{\\gtrsim}"),
                 ("≲", r"\\ensuremath{\\lesssim}"), ("ℰ", r"\\ensuremath{\\mathcal{E}}"),
                 ("ᵏ", r"\\ensuremath{^{k}}"), ("→", r"\\ensuremath{\\rightarrow}"),
                 ("√", r"\\ensuremath{\\surd}"),
                 ("⟨", r"\\ensuremath{\\langle}"), ("⟩", r"\\ensuremath{\\rangle}")]:
        p = re.sub(a, b, p)
    parts[i] = p
with tempfile.NamedTemporaryFile("w", suffix=".md", delete=False) as f:
    f.write("".join(parts)); tmp = f.name
out = os.path.join(HERE, "CMAGIC_COOKBOOK.pdf")
subprocess.run(["pandoc", tmp, "-o", out, "--pdf-engine=xelatex",
                "-V", "geometry:margin=25mm", "-V", "fontsize=11pt",
                "-V", "mainfont=STIX Two Text", "-V", "mathfont=STIX Two Math",
                "-V", "monofont=Menlo",
                "-V", "colorlinks=true", "--toc", "--toc-depth=2",
                "-M", "title=The CMAGIC Cookbook", "-M", "author=Lifan Wang"], check=True)
os.unlink(tmp)
print("wrote", out)
