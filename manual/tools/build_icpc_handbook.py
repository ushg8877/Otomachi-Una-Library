#!/usr/bin/env python3
"""Build generated-handbook.tex and generated/code; preserve the contest snapshot."""
from __future__ import annotations

import hashlib
import re
import pathlib

MANUAL_DIR = pathlib.Path(__file__).resolve().parents[1]
ROOT = MANUAL_DIR.parent
MANUAL = MANUAL_DIR / "manual.tex"
OUT_TEX = MANUAL_DIR / "generated-handbook.tex"
CODE_DIR = MANUAL_DIR / "generated" / "code"

BANNER = re.compile(
    r"^//{4,}\s*\r?\n(?://[^\n]*\r?\n)*?^//{4,}\s*\r?\n", re.MULTILINE
)
SEP = re.compile(
    r"^//\s*[-=\u2010-\u2015]{2,}[-=\s]*\s*\r?\n(?://[^\n]*\r?\n)*?"
    r"^//\s*[-=\u2010-\u2015]{2,}[-=\s]*\s*\r?\n",
    re.MULTILINE,
)
FOOTER = re.compile(r"\r?\n// end for[^\n]*\r?\n//+\s*$")

TEMPLATE_KEYS = (
    "template for",
    "usage:",
    "verify:",
    "fast mmap",
    "string algorithm template",
)

# manual order: (relative path, section label for LaTeX)
TEMPLATES: list[tuple[str, str]] = [
    ('basic/coordinate-compression.cpp', 'basic'),
    ('basic/fast-io.cpp', 'basic'),
    ('basic/template.cpp', 'basic'),
    ('data-structure/disjoint-set-rollback.cpp', 'data-structure'),
    ('data-structure/disjoint-set.cpp', 'data-structure'),
    ('data-structure/fast-set.cpp', 'data-structure'),
    ('data-structure/fenwick/fenwick.cpp', 'data-structure'),
    ('data-structure/fenwick/max.cpp', 'data-structure'),
    ('data-structure/fenwick/min.cpp', 'data-structure'),
    ('data-structure/fenwick/sum.cpp', 'data-structure'),
    ('data-structure/fhq-treap.cpp', 'data-structure'),
    ('data-structure/li-chao-tree.cpp', 'data-structure'),
    ('data-structure/li-chao-tree-static.cpp', 'data-structure'),
    ('data-structure/linear-rmq.cpp', 'data-structure'),
    ('data-structure/ordered-disjoint-interval-tree-fast.cpp', 'data-structure'),
    ('data-structure/ordered-disjoint-interval-tree.cpp', 'data-structure'),
    ('data-structure/priority-queue-two-stacks.cpp', 'data-structure'),
    ('data-structure/rollback-int.cpp', 'data-structure'),
    ('data-structure/rollback-int64.cpp', 'data-structure'),
    ('data-structure/segment-tree/dynamic.cpp', 'data-structure'),
    ('data-structure/segment-tree/iterative-prefix.cpp', 'data-structure'),
    ('data-structure/segment-tree/lazy.cpp', 'data-structure'),
    ('data-structure/segment-tree/merge.cpp', 'data-structure'),
    ('data-structure/segment-tree/point.cpp', 'data-structure'),
    ('data-structure/segment-tree/range-add-max.cpp', 'data-structure'),
    ('data-structure/segment-tree/range-add-min.cpp', 'data-structure'),
    ('data-structure/segment-tree/range-add-sum.cpp', 'data-structure'),
    ('data-structure/sparse-table.cpp', 'data-structure'),
    ('graph/block-cut-tree.cpp', 'graph'),
    ('graph/chordal.cpp', 'graph'),
    ('graph/dijkstra.cpp', 'graph'),
    ('graph/directed-mst.cpp', 'graph'),
    ('graph/edge-biconnected.cpp', 'graph'),
    ('graph/flow/max-flow-hlpp.cpp', 'graph'),
    ('graph/flow/min-cost-flow-isap.cpp', 'graph'),
    ('graph/flow/min-cost-flow-simplex.cpp', 'graph'),
    ('graph/matching/bipartite.cpp', 'graph'),
    ('graph/matching/general.cpp', 'graph'),
    ('graph/matching/regular-bipartite.cpp', 'graph'),
    ('graph/steiner-tree.cpp', 'graph'),
    ('graph/two-sat.cpp', 'graph'),
    ('tree/cartesian-tree.cpp', 'tree'),
    ('tree/centroid-decomposition.cpp', 'tree'),
    ('tree/divide-combine-tree.cpp', 'tree'),
    ('tree/heavy-light-decomposition.cpp', 'tree'),
    ('string/aho-corasick.cpp', 'string'),
    ('string/generalized-suffix-automaton.cpp', 'string'),
    ('string/kmp-z.cpp', 'string'),
    ('string/lyndon.cpp', 'string'),
    ('string/manacher.cpp', 'string'),
    ('string/palindromic-automaton.cpp', 'string'),
    ('string/rolling-hash.cpp', 'string'),
    ('string/runs.cpp', 'string'),
    ('string/suffix-array.cpp', 'string'),
    ('string/suffix-automaton.cpp', 'string'),
    ('math/big-int.cpp', 'math'),
    ('math/dynamic-mod-int.cpp', 'math'),
    ('math/gaussian-elimination.cpp', 'math'),
    ('math/linear-basis/field-prefix.cpp', 'math'),
    ('math/linear-basis/field.cpp', 'math'),
    ('math/linear-basis/xor-prefix.cpp', 'math'),
    ('math/linear-basis/xor.cpp', 'math'),
    ('math/mod-int.cpp', 'math'),
    ('math/rational.cpp', 'math'),
    ('math/xor-equations.cpp', 'math'),
    ('number-theory/dujiao-sieve.cpp', 'number-theory'),
    ('number-theory/floor-sum.cpp', 'number-theory'),
    ('number-theory/linear-sieve.cpp', 'number-theory'),
    ('number-theory/modular-arithmetic.cpp', 'number-theory'),
    ('number-theory/pollard-rho.cpp', 'number-theory'),
    ('number-theory/stern-brocot.cpp', 'number-theory'),
    ('polynomial/fwt.cpp', 'polynomial'),
    ('polynomial/integer.cpp', 'polynomial'),
    ('polynomial/mtt.cpp', 'polynomial'),
    ('polynomial/ntt-fast.cpp', 'polynomial'),
    ('polynomial/ntt-int64.cpp', 'polynomial'),
    ('polynomial/ntt.cpp', 'polynomial'),
    ('geometry/geometry-2d.cpp', 'geometry'),
]


def strip_leading_comment_run(text: str) -> str:
    lines = text.splitlines(keepends=True)
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if not s or s.startswith("//"):
            i += 1
            continue
        break
    return "".join(lines[i:])


def strip_template_source(text: str) -> str:
    text = FOOTER.sub("", text.rstrip())
    changed = True
    while changed:
        changed = False
        for pat in (BANNER, SEP):
            m = pat.search(text)
            while m:
                chunk = m.group(0)
                if any(k in chunk.lower() for k in TEMPLATE_KEYS):
                    text = text[: m.start()] + text[m.end() :]
                    changed = True
                    m = pat.search(text)
                else:
                    break
    text = strip_leading_comment_run(text)
    return text.strip() + "\n"


def prepare_doc_latex(s: str) -> str:
    return s.strip()


def braced_content(text: str, start: int) -> tuple[str, int]:
    assert text[start] == "{"
    depth = 0
    i = start
    while i < len(text):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[start + 1 : i], i + 1
        i += 1
    return "", len(text)


def parse_manual_docs(path: pathlib.Path) -> dict[str, dict[str, str]]:
    text = path.read_text(encoding="utf-8")
    docs: dict[str, dict[str, str]] = {}
    for m in re.finditer(r"\\subsection\{\\texttt\{", text):
        pos = m.end()
        name, pos = braced_content(text, pos - 1)
        name = name.replace("\\_", "_")
        block_end = text.find("\\subsection{", pos)
        if block_end < 0:
            block_end = text.find("\\section{", pos)
        if block_end < 0:
            block_end = len(text)
        block = text[pos:block_end]
        entry = {"usage": "", "api": "", "note": ""}
        for macro, key in (
            ("\\TemplUsage", "usage"),
            ("\\TemplApi", "api"),
            ("\\TemplNote", "note"),
        ):
            idx = block.find(macro)
            if idx < 0:
                continue
            j = idx + len(macro)
            if j < len(block) and block[j] == "[":
                j = block.find("]", j) + 1
            if j >= len(block) or block[j] != "{":
                continue
            content, _ = braced_content(block, j)
            entry[key] = prepare_doc_latex(content)
        # map filename key to relative paths
        keys = [name]
        if name.startswith("Segment Tree/"):
            keys.append(f"ds/{name}")
        elif "/" not in name and not name.endswith(".cpp"):
            pass
        else:
            for rel, _ in TEMPLATES:
                if rel.endswith("/" + name) or rel == name:
                    keys.append(rel)
        for k in keys:
            docs[k] = entry
            docs[k.split("/")[-1]] = entry
    return docs


def code_filename(rel: str) -> str:
    s = (
        rel.replace("/", "__")
        .replace(" ", "_")
        .replace("(", "_")
        .replace(")", "_")
    )
    if not s.isascii():
        s = "f_" + hashlib.sha1(rel.encode("utf-8")).hexdigest()[:12]
    if not s.endswith(".cpp"):
        s += ".cpp"
    return s


def latex_path(p: str) -> str:
    return p.replace("_", "\\_")


def write_handbook(docs: dict[str, dict[str, str]]) -> None:
    CODE_DIR.mkdir(parents=True, exist_ok=True)
    lines = [
        r"\documentclass[UTF8,a4paper,10pt]{ctexart}",
        r"\usepackage{amsmath,amssymb,geometry,listings,xcolor,hyperref,enumitem}",
        r"\geometry{margin=1.8cm}",
        r"\lstset{",
        r"  language=C++,",
        r"  basicstyle=\ttfamily\small,",
        r"  keywordstyle=\color{blue},",
        r"  commentstyle=\color{gray!70},",
        r"  frame=single,",
        r"  breaklines=true,",
        r"  tabsize=2,",
        r"  extendedchars=true,",
        r"  inputencoding=utf8,",
        r"  literate={—}{{-}}1,",
        r"}",
        r"\title{Otomachi Una ICPC 模板手册}",
        r"\author{Junlin Ye}",
        r"\date{\today}",
        r"\begin{document}",
        r"\maketitle",
        r"\tableofcontents",
        r"\newpage",
        r"\begin{center}\small 各节含中文使用说明与可粘贴代码（已去除文件首尾模板注释）。\end{center}",
        r"\bigskip",
    ]
    sec_labels = {
        "basic": "basic/",
        "ds": "ds/",
        "graph": "graph/",
        "trees": "trees/",
        "number theory": "number theory/",
        "polynomial": "polynomial/",
        "string": "string/",
        "geometry": "geometry/",
        "root": "根目录",
    }
    current_sec = None
    count = 0
    for rel, sec in TEMPLATES:
        src = ROOT / rel
        if not src.exists():
            print("MISSING", rel)
            continue
        code = strip_template_source(src.read_text(encoding="utf-8"))
        code_path = CODE_DIR / code_filename(rel)
        code_path.write_text(code, encoding="utf-8", newline="\n")
        if sec != current_sec:
            current_sec = sec
            lines.append(f"\\section{{\\texttt{{{latex_path(sec+'/')}}}}}")
        doc = docs.get(rel) or docs.get(rel.split("/")[-1], {"usage": "", "api": "", "note": ""})
        lines.append(f"\\subsection{{\\texttt{{{latex_path(rel)}}}}}")
        if doc["usage"]:
            lines.append("\\noindent\\textbf{使用说明：}" + doc["usage"] + "\\par")
        if doc["api"]:
            lines.append("\\noindent\\textbf{接口：}")
            lines.append(doc["api"])
            lines.append("\\par")
        if doc["note"]:
            lines.append("\\noindent\\textbf{注意事项：}" + doc["note"] + "\\par")
        lines.append("\\vspace{0.35em}")
        lst = str(code_path.relative_to(MANUAL_DIR)).replace("\\", "/")
        lines.append(f"\\lstinputlisting{{{lst}}}")
        lines.append("\\vspace{0.9em}")
        count += 1
    lines.extend([r"\end{document}", ""])
    OUT_TEX.write_text("\n".join(lines), encoding="utf-8", newline="\n")
    print(f"Wrote {OUT_TEX} ({count} templates)")


if __name__ == "__main__":
    write_handbook(parse_manual_docs(MANUAL))
