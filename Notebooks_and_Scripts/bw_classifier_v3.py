"""
Ballesteros-Weinstein (BW) Classifier for GPCR Variants  v3.0
==============================================================
Supports TWO residue table formats from GPCRdb:

  FORMAT 1 (residue_table.xlsx  — original 19 receptors):
    Columns:  BW  |  Receptor1  |  Receptor2  | ...
    BW cell:  "1x20 1.20"  (GPCRdb + BW notation combined)
    Segment headers appear in the BW column.

  FORMAT 2 (residue_table2.xlsx — new receptors):
    Columns:  GPCRdb(A)  |  BW  |  Receptor1  |  Receptor2  | ...
    col0:     "1x20"     col1: "1.20"
    Segment headers appear in col0 only.

The script auto-detects which format is present and parses accordingly.
Both files can be passed together via --residue flag (comma-separated).

USAGE
-----
  Single residue table:
    python bw_classifier_v3.py <variants_folder> <residue_table.xlsx> <output_folder>

  Multiple residue tables (combined):
    python bw_classifier_v3.py <variants_folder> "table1.xlsx,table2.xlsx" <output_folder>

VARIANT CSV requirements
------------------------
  Must contain columns (case-insensitive, stripped):
    protein position   Amino acid change   consequences

  Filename must contain the receptor name or a known alias, e.g.:
    adra1a_variants.csv        → ADRA1A
    adrb2_exon_coding.csv      → ADRB2
    ffar1_exon_coding.csv      → FFAR1

SUPPORTED RECEPTORS
-------------------
  From residue_table.xlsx  (19 receptors):
    FFAR1  FFAR2  FFAR3  GLP1R  GLP2R  GPR142 GPR183 GPR75
    GPBAR  LPAR1  LPAR3  LPAR5  LPAR6  MC4R   S1PR1  S1PR2
    S1PR5  SCTR   VIPR1

  From residue_table2.xlsx (12 new receptors):
    ADRA1A  ADRB2   ADRB3   AGTR1   CASR    CNR1
    GHSR    MTNR1B  NPY1R   PTGER3  SSTR2   SUCNR1
"""

import os, sys, re, glob
import pandas as pd
from openpyxl import load_workbook

# ── Receptor metadata ─────────────────────────────────────────────────────────
RECEPTOR_META = {
    # Original 19
    "FFAR1": {"class":"A",  "g_protein":"Gq/11, Gi",     "uniprot":"O14842"},
    "FFAR2": {"class":"A",  "g_protein":"Gi, Gq",        "uniprot":"O15552"},
    "FFAR3": {"class":"A",  "g_protein":"Gi",             "uniprot":"O14843"},
    "GLP1R": {"class":"B1", "g_protein":"Gs",             "uniprot":"P43220"},
    "GLP2R": {"class":"B1", "g_protein":"Gs",             "uniprot":"O95838"},
    "GPR142":{"class":"A",  "g_protein":"Gq/11",          "uniprot":"Q7Z601"},
    "GPR183":{"class":"A",  "g_protein":"Gi",             "uniprot":"P32249"},
    "GPR75": {"class":"A",  "g_protein":"Gq/11",          "uniprot":"Q9BZJ8"},
    "GPBAR": {"class":"A",  "g_protein":"Gs",             "uniprot":"Q8TDU6"},
    "LPAR1": {"class":"A",  "g_protein":"Gi, Gq, G12/13","uniprot":"Q92633"},
    "LPAR3": {"class":"A",  "g_protein":"Gi, Gq",         "uniprot":"Q9UBY5"},
    "LPAR5": {"class":"A",  "g_protein":"G12/13, Gq",     "uniprot":"Q9H1C0"},
    "LPAR6": {"class":"A",  "g_protein":"Gs, Gi",         "uniprot":"P43657"},
    "MC4R":  {"class":"A",  "g_protein":"Gs, Gq",         "uniprot":"P32245"},
    "S1PR1": {"class":"A",  "g_protein":"Gi",             "uniprot":"P21453"},
    "S1PR2": {"class":"A",  "g_protein":"G12/13, Gi, Gq","uniprot":"O95136"},
    "S1PR5": {"class":"A",  "g_protein":"Gi",             "uniprot":"Q9H228"},
    "SCTR":  {"class":"B1", "g_protein":"Gs",             "uniprot":"P47872"},
    "VIPR1": {"class":"B1", "g_protein":"Gs",             "uniprot":"P32241"},
    # New 12
    "ADRA1A":{"class":"A",  "g_protein":"Gq/11",          "uniprot":"P35348"},
    "ADRB2": {"class":"A",  "g_protein":"Gs",             "uniprot":"P07550"},
    "ADRB3": {"class":"A",  "g_protein":"Gs",             "uniprot":"P13945"},
    "AGTR1": {"class":"A",  "g_protein":"Gq/11, G12/13",  "uniprot":"P30556"},
    "CASR":  {"class":"C",  "g_protein":"Gq/11, Gi",      "uniprot":"P41180"},
    "CNR1":  {"class":"A",  "g_protein":"Gi",             "uniprot":"P21554"},
    "GHSR":  {"class":"A",  "g_protein":"Gq/11",          "uniprot":"Q92847"},
    "MTNR1B":{"class":"A",  "g_protein":"Gi",             "uniprot":"P49286"},
    "NPY1R": {"class":"A",  "g_protein":"Gi",             "uniprot":"P25929"},
    "PTGER3":{"class":"A",  "g_protein":"Gi",             "uniprot":"P43115"},
    "SSTR2": {"class":"A",  "g_protein":"Gi",             "uniprot":"P30874"},
    "SUCNR1":{"class":"A",  "g_protein":"Gi, Gq/11",      "uniprot":"Q9BXA5"},
    # Additional receptor found in Remaining_GPCRs folder
    "CHRM3": {"class":"A",  "g_protein":"Gq/11",          "uniprot":"P20309"},
}

# ── Functional notes for BW positions ────────────────────────────────────────
CLASS_A_NOTES = {
    "1.50":"N1.50 — structural pivot TM1; Na+/water channel",
    "2.50":"D2.50 — Na+ binding pocket; allosteric modulation",
    "3.49":"D/E3.49 — DRY/ERY motif; ionic lock partner",
    "3.50":"R3.50 — DRY motif; direct Gα contact upon activation",
    "3.51":"Y/H3.51 — DRY motif; ionic lock",
    "4.50":"W4.50 — structural pivot TM4",
    "5.50":"P5.50 — structural pivot TM5; flexibility for outward movement",
    "6.47":"C6.47 — CWxP motif",
    "6.48":"W6.48 — CWxP toggle switch; rotates upon activation",
    "6.50":"P6.50 — CWxP motif; TM6 kink",
    "7.49":"N7.49 — NPxxY motif",
    "7.50":"P7.50 — NPxxY motif; structural pivot TM7",
    "7.53":"Y7.53 — NPxxY motif; reorients toward TM3 on activation",
    "3.44":"Upper TM3 — orthosteric ligand binding pocket",
    "5.41":"TM5 — ligand binding pocket",
    "5.46":"TM5 — activation microswitch (4.53–5.46 contact pair)",
    "6.33":"TM6 cytoplasmic — direct G-protein α5 helix contact",
    "7.60":"Lower TM7 — proximal to H8; Gα contact region",
}
CLASS_B_NOTES = {
    "2.46":"R2.46 — conserved class B Gαs contact",
    "3.50":"Structural pivot TM3",
    "3.53":"Y3.53 — conserved class B Gαs contact",
    "3.54":"L3.54 — conserved class B Gαs contact",
    "3.58":"L/I3.58 — conserved class B Gαs contact",
    "5.50":"Structural pivot TM5",
    "5.64":"K5.64 — direct H-bond to Gα α5 helix (Q384/R385 of Gαs)",
    "7.50":"Structural pivot TM7",
    "7.61":"N7.61 — TM7–H8 junction; Gαs contact",
}
RECEPTOR_NOTES = {
    "FFAR1": {"6.33":"DIRECT Gq contact — cryo-EM confirmed (PNAS 2023)",
              "6.48":"V6.48 (not W) — atypical CWxP toggle switch in FFAR1",
              "7.53":"T7.53 (not Y) — atypical NPxxT in FFAR1"},
    "GPBAR": {"3.49":"E3.49 — ERY variant (not DRY) in GPBAR/TGR5 (PDB 7CFN)"},
    "MC4R":  {"3.50":"R3.50 = R131 — DRY (PDB 7AUE)",
              "7.53":"Y7.53 = Y265 — NPxxY"},
    "ADRB2": {"3.50":"R3.50 — DRY motif (PDB 3SN6 β2AR-Gs)",
              "6.48":"W6.48 — CWxP toggle switch (PDB 3SN6)",
              "7.53":"Y7.53 — NPxxY motif"},
    "CNR1":  {"3.50":"R3.50 — DRY motif (PDB 5XRA CB1-Gi)",
              "6.48":"W6.48 — CWxP toggle switch"},
    "AGTR1": {"3.50":"R3.50 — DRY motif (PDB 6OS0 AT1R)",
              "6.48":"W6.48 — CWxP toggle switch"},
    "GLP1R": {"5.64":"K5.64 = K334 — H-bond to Gα α5 (PDB 6X18)"},
    "SCTR":  {"5.64":"K5.64 = K319 — H-bond to Gα α5 (PDB 7EIY)"},
    "VIPR1": {"5.64":"K5.64 = K330 — H-bond to Gα α5 (PDB 6VN7)"},
    "CHRM3": {"3.50":"R3.50 — DRY motif; critical for Gq/11 coupling",
              "6.48":"W6.48 — CWxP toggle switch"},
}

# ── Filename alias detection ──────────────────────────────────────────────────
ALIASES = {
    # original 19
    "ffar1":"FFAR1","gpr40":"FFAR1",
    "ffar2":"FFAR2","gpr43":"FFAR2",
    "ffar3":"FFAR3","gpr41":"FFAR3",
    "glp1r":"GLP1R","glp2r":"GLP2R",
    "gpr142":"GPR142","gpr183":"GPR183","ebi2":"GPR183",
    "gpr75":"GPR75","gpbar":"GPBAR","tgr5":"GPBAR","gpbar1":"GPBAR",
    "lpar1":"LPAR1","lpar3":"LPAR3","lpar5":"LPAR5","lpar6":"LPAR6",
    "mc4r":"MC4R",
    "s1pr1":"S1PR1","edg1":"S1PR1","s1p1":"S1PR1",
    "s1pr2":"S1PR2","edg5":"S1PR2","s1p2":"S1PR2",
    "s1pr5":"S1PR5","edg8":"S1PR5","s1p5":"S1PR5",
    "sctr":"SCTR","secretin":"SCTR",
    "vipr1":"VIPR1","vpac1":"VIPR1","vip1r":"VIPR1",
    # new 12
    "adra1a":"ADRA1A","adra1":"ADRA1A","alpha1a":"ADRA1A",
    "adrb2":"ADRB2","b2ar":"ADRB2","beta2":"ADRB2",
    "adrb3":"ADRB3","b3ar":"ADRB3","beta3":"ADRB3",
    "agtr1":"AGTR1","at1r":"AGTR1","at1":"AGTR1",
    "casr":"CASR","cas":"CASR","calcium":"CASR",
    "cnr1":"CNR1","cb1":"CNR1","cb1r":"CNR1",
    "ghsr":"GHSR","ghrelin":"GHSR",
    "mtnr1b":"MTNR1B","mt2":"MTNR1B","mt2r":"MTNR1B","melatonin":"MTNR1B",
    "npy1r":"NPY1R","y1":"NPY1R","y1r":"NPY1R",
    "ptger3":"PTGER3","ep3":"PTGER3","ep3r":"PTGER3",
    "sstr2":"SSTR2","sst2":"SSTR2","somatostatin":"SSTR2",
    "sucnr1":"SUCNR1","gpr91":"SUCNR1","succinate":"SUCNR1",
    "chrm3":"CHRM3","m3":"CHRM3","muscarinic3":"CHRM3","m3r":"CHRM3",
}

def detect_receptor(filename):
    base = os.path.basename(filename).lower()

    # ── Step 1: stem-based detection (most reliable) ──────────────────────
    # Handles both:
    #   ADRA1A_coding_variants.csv         (new Remaining_GPCRs pattern)
    #   ffar1_exon_coding_variants.csv     (original pattern)
    m = re.match(r'^(.+?)(?:_exon)?_coding_variants\.csv$', base)
    if m:
        stem = m.group(1).upper()
        if stem in RECEPTOR_META:          # exact canonical name match
            return stem
        # stem might itself be an alias (e.g. "gpr40" for FFAR1)
        if stem.lower() in ALIASES:
            return ALIASES[stem.lower()]

    # ── Step 2: alias substring search (fallback for other naming styles) ──
    for alias in sorted(ALIASES, key=len, reverse=True):
        if alias in base:
            return ALIASES[alias]

    return None

# ── Segment headers ───────────────────────────────────────────────────────────
SEGMENT_HEADERS = {
    "TM1","ICL1","TM2","ECL1","TM3","ICL2",
    "TM4","ECL2","TM5","ICL3","TM6","ECL3","TM7","H8"
}
SEGMENT_ORDER = ["TM1","ICL1","TM2","ECL1","TM3","ICL2",
                 "TM4","ECL2","TM5","ICL3","TM6","ECL3","TM7","H8"]

# ── Format detection ──────────────────────────────────────────────────────────
def detect_format(rows):
    """
    Format 1: header[0] contains 'BW' or similar, residue columns start at index 1
    Format 2: header[0]='GPCRdb(A)', header[1]='BW', residue columns start at index 2
    """
    h = rows[0]
    h0 = str(h[0]).strip().lower() if h[0] else ""
    h1 = str(h[1]).strip().lower() if len(h)>1 and h[1] else ""
    if "gpcrdb" in h0 or ("bw" in h1 and "bw" not in h0):
        return 2   # new format
    return 1       # original format

# ── Residue table parser ──────────────────────────────────────────────────────
def parse_residue_table(xlsx_path):
    """
    Parse either format and return:
      pos_map[receptor][seq_pos] = {"bw": "1.50", "segment": "TM1"}
      seg_bounds[receptor][segment] = (min_pos, max_pos)
    Also returns col_to_rec dict for this file.
    """
    wb = load_workbook(xlsx_path, read_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))
    fmt  = detect_format(rows)

    header      = rows[0]
    rec_col_start = 2 if fmt == 2 else 1
    bw_col        = 1 if fmt == 2 else 0   # column that holds BW number

    # Map column index → canonical receptor name
    # For format 2, header names are long like "α1A-adrenoceptor Human"
    # Build a reverse lookup from display name → canonical
    DISPLAY_TO_CANONICAL = {
        "α1A-adrenoceptor Human": "ADRA1A",
        "β2-adrenoceptor Human":  "ADRB2",
        "β3-adrenoceptor Human":  "ADRB3",
        "AT1 receptor Human":     "AGTR1",
        "CaS receptor Human":     "CASR",
        "CB1 receptor Human":     "CNR1",
        "ghrelin receptor Human": "GHSR",
        "MT2 receptor Human":     "MTNR1B",
        "Y1 receptor Human":      "NPY1R",
        "EP3 receptor Human":     "PTGER3",
        "SST2 receptor Human":    "SSTR2",
        "succinate receptor Human":"SUCNR1",
        # Format 1 headers (short names already in ALIASES)
        "FFA1 receptor Human":    "FFAR1",
        "FFA2 receptor Human":    "FFAR2",
        "FFA3 receptor Human":    "FFAR3",
        "GLP-1 receptor Human":   "GLP1R",
        "GLP-2 receptor Human":   "GLP2R",
        "GPR142 Human":           "GPR142",
        "GPR183 Human":           "GPR183",
        "GPR75 Human":            "GPR75",
        "GPBA receptor Human":    "GPBAR",
        "LPA1 receptor Human":    "LPAR1",
        "LPA3 receptor Human":    "LPAR3",
        "LPA5 receptor Human":    "LPAR5",
        "LPA6 receptor Human":    "LPAR6",
        "MC4 receptor Human":     "MC4R",
        "S1P1 receptor Human":    "S1PR1",
        "S1P2 receptor Human":    "S1PR2",
        "S1P5 receptor Human":    "S1PR5",
        "secretin receptor Human":"SCTR",
        "VPAC1 receptor Human":   "VIPR1",
    }

    col_to_rec = {}
    for ci in range(rec_col_start, len(header)):
        h = str(header[ci]).strip() if header[ci] else ""
        if h in DISPLAY_TO_CANONICAL:
            col_to_rec[ci] = DISPLAY_TO_CANONICAL[h]
        elif h:
            # Try alias detection on the display name
            rec = detect_receptor(h.lower().replace(" ","_"))
            if rec:
                col_to_rec[ci] = rec

    pos_map   = {rec: {} for rec in col_to_rec.values()}
    seg_pos   = {rec: {} for rec in col_to_rec.values()}
    current_seg = "N-term"

    for row in rows[1:]:
        # Segment header detection: look in col0 for format2, bw_col for format1
        seg_check = str(row[0]).strip() if row[0] else ""
        if seg_check in SEGMENT_HEADERS:
            current_seg = seg_check
            continue

        # Extract BW number
        bw_raw = row[bw_col]
        if not bw_raw:
            continue
        bw_str = str(bw_raw).strip()
        # Handle combined cells like "1x20 1.20" or plain "1.20"
        bw_matches = re.findall(r'\d+\.\d+', bw_str)
        if not bw_matches:
            continue
        bw_num = bw_matches[0]

        # Parse each receptor column
        for ci, rec_name in col_to_rec.items():
            cell_val = row[ci]
            if not cell_val or str(cell_val).strip() in ('-','','None'):
                continue
            m = re.match(r'^([A-Z]+)(\d+)$', str(cell_val).strip())
            if not m:
                continue
            seq_pos = int(m.group(2))
            pos_map[rec_name][seq_pos] = {
                "bw": bw_num, "segment": current_seg, "aa": m.group(1)
            }
            seg_pos[rec_name].setdefault(current_seg, []).append(seq_pos)

    # Segment bounds
    seg_bounds = {}
    for rec, segs in seg_pos.items():
        seg_bounds[rec] = {seg: (min(p), max(p)) for seg, p in segs.items()}

    print(f"  [{os.path.basename(xlsx_path)}] Format {fmt} → "
          f"{len(col_to_rec)} receptors, "
          f"{sum(len(v) for v in pos_map.values())} residues")
    return pos_map, seg_bounds

# ── Structural class inference ────────────────────────────────────────────────
def infer_segment(pos, seg_bounds):
    tm_segs = {s: b for s, b in seg_bounds.items() if s in SEGMENT_HEADERS}
    if not tm_segs:
        return "Unknown"
    sorted_segs = sorted(tm_segs.items(), key=lambda x: x[1][0])
    first_start = sorted_segs[0][1][0]
    last_end    = sorted_segs[-1][1][1]
    if pos < first_start: return "N-term"
    if pos > last_end:    return "C-term"
    for i in range(len(sorted_segs)-1):
        seg_a, (_, end_a) = sorted_segs[i]
        seg_b, (start_b, _) = sorted_segs[i+1]
        if end_a < pos < start_b:
            idx_a = SEGMENT_ORDER.index(seg_a) if seg_a in SEGMENT_ORDER else -1
            if 0 <= idx_a+1 < len(SEGMENT_ORDER):
                return SEGMENT_ORDER[idx_a+1]
            return f"between {seg_a}-{seg_b}"
    return "Unknown"

def structural_class(segment):
    if segment.startswith("TM"):  return "Transmembrane"
    if segment.startswith("ICL"): return "Intracellular loop"
    if segment.startswith("ECL"): return "Extracellular loop"
    if segment == "H8":           return "Helix 8"
    if segment == "N-term":       return "N-terminus"
    if segment == "C-term":       return "C-terminus"
    return segment

def get_note(bw, receptor, rec_class):
    if bw == "—": return ""
    rec_notes = RECEPTOR_NOTES.get(receptor, {})
    if bw in rec_notes: return rec_notes[bw]
    pool = CLASS_B_NOTES if rec_class in ("B1","B") else CLASS_A_NOTES
    return pool.get(bw, "")

# ── Classify one variant file ─────────────────────────────────────────────────
# Column name variants to try for each required field
COL_CANDIDATES = {
    "protein position": [
        "protein position", "protein_position", "position",
        "protein pos", "pos", "aa position", "amino acid position",
        "protein.position", "prot_position", "seq position",
    ],
    "Amino acid change": [
        "amino acid change", "amino_acid_change", "aa change",
        "aa_change", "aa substitution", "change", "mutation",
        "variant", "residue change",
    ],
    "consequences": [
        "consequences", "consequence", "variant type",
        "variant_type", "effect", "variant effect",
    ],
}

def find_col(columns, field):
    """Return the actual column name matching the requested field (case-insensitive)."""
    cols_lower = {c.lower().strip(): c for c in columns}
    for candidate in COL_CANDIDATES.get(field, [field]):
        if candidate.lower() in cols_lower:
            return cols_lower[candidate.lower()]
    # Last resort: partial match
    for col_low, col_orig in cols_lower.items():
        if field.lower().replace(" ","") in col_low.replace(" ","").replace("_",""):
            return col_orig
    return None


def classify_file(variant_csv, pos_map, seg_bounds, receptor_name):
    meta      = RECEPTOR_META[receptor_name]
    rec_map   = pos_map.get(receptor_name, {})
    rec_bounds = seg_bounds.get(receptor_name, {})

    df = pd.read_csv(variant_csv)
    df.columns = [c.strip() for c in df.columns]

    # Detect actual column names in this file
    pos_col  = find_col(df.columns, "protein position")
    aa_col   = find_col(df.columns, "Amino acid change")
    cons_col = find_col(df.columns, "consequences")

    if pos_col is None:
        raise ValueError(
            f"Cannot find a 'protein position' column in {os.path.basename(variant_csv)}\n"
            f"  Columns found: {list(df.columns)}\n"
            f"  Expected one of: {COL_CANDIDATES['protein position']}"
        )

    out_rows = []
    for _, row in df.iterrows():
        r = row.copy()
        try:
            pos = int(row[pos_col])
        except (ValueError, TypeError):
            r["Receptor"]           = receptor_name
            r["Structural region"]  = "Unknown"
            r["BW number"]          = "—"
            r["Structural class"]   = "Unknown"
            r["Receptor class"]     = meta["class"]
            r["G-protein coupling"] = meta["g_protein"]
            r["UniProt"]            = meta["uniprot"]
            r["Functional note"]    = "Non-numeric position"
            out_rows.append(r); continue

        if pos in rec_map:
            entry   = rec_map[pos]
            bw      = entry["bw"]
            segment = entry["segment"]
        else:
            bw      = "—"
            segment = infer_segment(pos, rec_bounds)

        r["Receptor"]           = receptor_name
        r["Structural region"]  = segment
        r["BW number"]          = bw
        r["Structural class"]   = structural_class(segment)
        r["Receptor class"]     = meta["class"]
        r["G-protein coupling"] = meta["g_protein"]
        r["UniProt"]            = meta["uniprot"]
        r["Functional note"]    = get_note(bw, receptor_name, meta["class"])
        out_rows.append(r)

    return pd.DataFrame(out_rows)

# ── Batch runner ──────────────────────────────────────────────────────────────
def run(variants_folder, xlsx_paths_str, output_folder):
    xlsx_paths = [p.strip() for p in xlsx_paths_str.split(",")]

    print("Loading residue tables...")
    combined_pos_map    = {}
    combined_seg_bounds = {}
    for xpath in xlsx_paths:
        pm, sb = parse_residue_table(xpath)
        combined_pos_map.update(pm)
        combined_seg_bounds.update(sb)

    total_receptors = len([r for r,v in combined_pos_map.items() if v])
    total_residues  = sum(len(v) for v in combined_pos_map.values())
    print(f"  Total: {total_receptors} receptors, {total_residues} residues loaded.\n")

    os.makedirs(output_folder, exist_ok=True)
    csv_files = glob.glob(os.path.join(variants_folder, "*.csv"))
    if not csv_files:
        print(f"No CSV files found in: {variants_folder}"); return

    all_dfs, ok, failed = [], 0, []

    for vpath in sorted(csv_files):
        rec = detect_receptor(vpath)
        if rec is None:
            print(f"  [SKIP] Cannot detect receptor from: {os.path.basename(vpath)}")
            failed.append(vpath); continue
        if rec not in RECEPTOR_META:
            print(f"  [SKIP] {rec} not in metadata table.")
            failed.append(vpath); continue
        if not combined_pos_map.get(rec):
            print(f"  [SKIP] {rec} not found in any residue table.")
            failed.append(vpath); continue
        try:
            df  = classify_file(vpath, combined_pos_map, combined_seg_bounds, rec)
            out = os.path.join(output_folder, f"{rec}_BW_classified.csv")
            df.to_csv(out, index=False)
            n_bw   = (df["BW number"] != "—").sum()
            n_loop = (df["BW number"] == "—").sum()
            print(f"  [OK]  {rec:10s} — {len(df):3d} variants "
                  f"({n_bw} with BW#, {n_loop} in loops/termini)")
            all_dfs.append(df); ok += 1
        except Exception as e:
            import traceback
            print(f"  [ERR] {rec}: {e}")
            traceback.print_exc()
            failed.append(vpath)

    if all_dfs:
        combined = pd.concat(all_dfs, ignore_index=True)
        cpath    = os.path.join(output_folder, "ALL_receptors_BW_classified.csv")
        combined.to_csv(cpath, index=False)
        print(f"\n  Combined → {cpath}")
        print(f"  {len(combined)} total variants across {ok} receptors")

    if failed:
        print(f"\n  Skipped/failed: {[os.path.basename(f) for f in failed]}")

# ── Entry point ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) == 4:
        run(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print(__doc__)
