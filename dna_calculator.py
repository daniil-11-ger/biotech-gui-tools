import streamlit as st
import pandas as pd
import numpy as np
from io import StringIO
from collections import Counter
import json
import urllib.parse

# --- CONSTANTS AND DATA ---
GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

AA_WEIGHTS = {
    'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2, 'E': 147.1, 'Q': 146.1, 
    'G': 75.1, 'H': 155.2, 'I': 131.2, 'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 
    'P': 115.1, 'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1, '_': 0
}

AA_pKa = {
    'D': 3.9, 'E': 4.3, 'C': 8.3, 'Y': 10.1, 'H': 6.0,
    'K': 10.5, 'R': 12.5, 'N_term': 8.0, 'C_term': 3.1
}

HYDROPATHY_INDEX = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 
    'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

RESTRICTION_ENZYMES = {
    'EcoRI': ('GAATTC', 'G^AATTC'),
    'HindIII': ('AAGCTT', 'A^AGCTT'),
    'BamHI': ('GGATCC', 'G^GATCC'),
    'NotI': ('GCGGCCGC', 'GC^GGCCGC'),
    'XhoI': ('CTCGAG', 'C^TCGAG'),
    'SalI': ('GTCGAC', 'G^TCGAC'),
    'PstI': ('CTGCAG', 'CTGCA^G'),
    'SmaI': ('CCCGGG', 'CCC^GGG'),
    'KpnI': ('GGTACC', 'GGTAC^C'),
    'SacI': ('GAGCTC', 'GAGCT^C')
}

# --- ANALYSIS FUNCTIONS ---
def get_reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(dna))

def calculate_protein_mass(protein):
    return sum(AA_WEIGHTS.get(aa, 0) for aa in protein)

def calculate_charge_at_ph_simple(protein_seq, ph):
    """Simple protein charge calculation"""
    if not protein_seq:
        return 0
    
    charge = 0.0
    
    # N-terminus
    if ph < 8.0:
        charge += 1
    
    # C-terminus
    if ph > 3.0:
        charge -= 1
    
    # Side chains
    for aa in protein_seq:
        if aa in ['D', 'E']:  # Acidic
            if ph > 4.0:
                charge -= 1
        elif aa in ['K', 'R']:  # Basic
            if ph < 10.0:
                charge += 1
        elif aa == 'H':  # Histidine
            if ph < 6.5:
                charge += 1
    
    return charge

def calculate_isoelectric_point_simple(protein_seq):
    """Simple pI calculation"""
    if not protein_seq:
        return 7.0
    
    acidic = protein_seq.count('D') + protein_seq.count('E')
    basic = protein_seq.count('R') + protein_seq.count('K') + protein_seq.count('H')
    
    if acidic + basic == 0:
        return 6.5
    
    # Simple formula
    pI = 6.5 + (basic - acidic) * 0.15
    return max(3.0, min(11.0, pI))

def find_restriction_sites(dna, enzymes=None):
    """Find restriction sites"""
    if enzymes is None:
        enzymes = RESTRICTION_ENZYMES
    
    sites = []
    for name, (seq, cut_site) in enzymes.items():
        pos = 0
        while pos < len(dna):
            found_pos = dna.find(seq, pos)
            if found_pos == -1:
                break
            sites.append({
                "Enzyme": name,
                "Sequence": seq,
                "Position": found_pos + 1,
                "Cut Site": cut_site,
                "End": found_pos + len(seq)
            })
            pos = found_pos + 1
    
    return sorted(sites, key=lambda x: x["Position"])

def calculate_hydropathy_profile(protein_seq, window=9):
    """Calculate hydropathy profile"""
    if len(protein_seq) < window:
        return [], []
    
    positions = []
    values = []
    
    for i in range(len(protein_seq) - window + 1):
        window_seq = protein_seq[i:i+window]
        hydropathy = np.mean([HYDROPATHY_INDEX.get(aa, 0) for aa in window_seq])
        positions.append(i + window//2)
        values.append(hydropathy)
    
    return positions, values

def find_orfs(dna, min_length=30):
    """Find open reading frames"""
    orfs = []
    rev_dna = get_reverse_complement(dna)
    
    for strand_name, strand in [('+', dna), ('-', rev_dna)]:
        for frame in range(3):
            seq = strand[frame:]
            for i in range(0, len(seq)-2, 3):
                codon = seq[i:i+3]
                if codon == 'ATG':  # Start codon
                    protein = []
                    for j in range(i, len(seq)-2, 3):
                        codon = seq[j:j+3]
                        aa = GENETIC_CODE.get(codon, 'X')
                        if aa == '_':  # Stop codon
                            protein.append('_')
                            orf_length = j + 3 - i
                            if orf_length >= min_length:
                                orfs.append({
                                    'Strand': strand_name,
                                    'Frame': frame + 1,
                                    'Start': i + 1,
                                    'End': j + 3,
                                    'Length': orf_length,
                                    'Protein': ''.join(protein)
                                })
                            break
                        protein.append(aa)
    
    return sorted(orfs, key=lambda x: x['Length'], reverse=True)

def calculate_gc_skew(dna, window=100):
    """Calculate GC-skew"""
    if len(dna) < window:
        return [], []
    
    positions = []
    skews = []
    
    for i in range(0, len(dna) - window + 1, window//2):
        window_seq = dna[i:i+window]
        g_count = window_seq.count('G')
        c_count = window_seq.count('C')
        total = g_count + c_count
        
        if total > 0:
            skew = (g_count - c_count) / total
        else:
            skew = 0
        
        positions.append(i + window//2)
        skews.append(skew)
    
    return positions, skews

def create_text_plasmid_map(dna_seq, sites):
    """Create text plasmid map"""
    if not sites:
        return "No restriction sites found"
    
    map_text = f"Plasmid Map ({len(dna_seq)} bp)\n"
    map_text += "0" + "-" * 50 + f"{len(dna_seq)} bp\n"
    
    for site in sorted(sites, key=lambda x: x['Position']):
        pos = site['Position']
        map_pos = int((pos / len(dna_seq)) * 50)
        map_text += " " * map_pos + "^\n"
        map_text += " " * map_pos + f"{site['Enzyme']} ({pos} bp)\n"
    
    return map_text

def create_simple_charge_plot(protein_seq):
    """Create simple charge plot"""
    pH_values = list(range(0, 15))
    charges = [calculate_charge_at_ph_simple(protein_seq, pH) for pH in pH_values]
    
    # Create text plot
    plot_text = "Charge vs pH:\n"
    plot_text += "pH: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14\n"
    plot_text += "Charge: "
    
    for charge in charges:
        if charge > 0:
            plot_text += f"+{int(charge)} "
        else:
            plot_text += f"{int(charge)} "
    
    # Find pI (where charge is closest to 0)
    pI = None
    min_charge = float('inf')
    for pH, charge in zip(pH_values, charges):
        if abs(charge) < abs(min_charge):
            min_charge = charge
            pI = pH
    
    plot_text += f"\n\nIsoelectric point (pI) ≈ {pI}"
    return plot_text

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    """Needleman-Wunsch alignment algorithm"""
    n, m = len(seq1), len(seq2)
    
    # Initialize matrix
    score = np.zeros((n+1, m+1))
    for i in range(1, n+1):
        score[i][0] = i * gap
    for j in range(1, m+1):
        score[0][j] = j * gap
    
    # Fill matrix
    for i in range(1, n+1):
        for j in range(1, m+1):
            match_score = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete = score[i-1][j] + gap
            insert = score[i][j-1] + gap
            score[i][j] = max(match_score, delete, insert)
    
    # Traceback
    align1, align2 = "", ""
    i, j = n, m
    
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score[i][j] == score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif i > 0 and score[i][j] == score[i-1][j] + gap:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1
    
    # Statistics
    matches = sum(1 for a, b in zip(align1, align2) if a == b and a != '-')
    identity = matches / max(len(seq1), len(seq2)) * 100
    
    return {
        'score': score[n][m],
        'identity': identity,
        'alignment1': align1,
        'alignment2': align2,
        'matches': matches,
        'gaps': align1.count('-') + align2.count('-')
    }

# --- INTERFACE ---
st.set_page_config(page_title="BioTech OS Pro", layout="wide")

# Custom CSS
st.markdown("""
    <style>
    .stApp { 
        background: linear-gradient(135deg, #0B0E14 0%, #161B22 100%);
        color: #FFFFFF !important; 
    }
    
    h1, h2, h3, h4 { 
        color: #10B981 !important; 
        font-family: 'Segoe UI', sans-serif;
    }
    
    .stTextArea textarea { 
        background-color: #1E2229 !important; 
        color: #00FF88 !important; 
        font-family: 'JetBrains Mono', monospace;
        border: 1px solid #30363D;
        border-radius: 8px;
    }
    
    .metric-box {
        background: linear-gradient(145deg, #161B22 0%, #1E2229 100%);
        border: 1px solid #30363D;
        border-radius: 12px;
        padding: 15px;
        text-align: center;
        margin-bottom: 10px;
    }
    
    .metric-val { 
        color: #10B981; 
        font-size: 1.5rem; 
        font-weight: bold;
    }
    
    .stTabs [data-baseweb="tab"] { 
        background-color: #1E2229;
        color: #8B949E !important;
        border-radius: 8px 8px 0 0;
        padding: 10px 20px;
    }
    
    .stTabs [aria-selected="true"] { 
        background-color: #161B22 !important;
        color: #10B981 !important;
        border-bottom: 3px solid #10B981 !important;
    }
    
    .stButton > button, .stDownloadButton > button {
        color: white !important;
        background-color: #10B981 !important;
        border: none !important;
        font-weight: bold !important;
        border-radius: 8px !important;
    }
    
    .stButton > button:hover, .stDownloadButton > button:hover {
        background-color: #0D8E6A !important;
    }
    
    .info-box {
        background-color: #1E2229;
        border-left: 4px solid #10B981;
        padding: 15px;
        margin: 10px 0;
        border-radius: 0 8px 8px 0;
    }
    
    .alignment-match { color: #10B981; font-weight: bold; }
    .alignment-mismatch { color: #FF6B6B; }
    .alignment-gap { color: #8B949E; }
    
    .blast-container {
        background-color: #1E2229;
        padding: 15px;
        border-radius: 10px;
        margin: 10px 0;
    }
    
    .export-section {
        background-color: #1E2229;
        padding: 20px;
        border-radius: 10px;
        margin: 15px 0;
    }
    
    .simple-plot {
        font-family: 'Courier New', monospace;
        background-color: #0B0E14;
        color: #00FF88;
        padding: 15px;
        border-radius: 8px;
        overflow-x: auto;
        white-space: pre;
        font-size: 12px;
    }
    </style>
""", unsafe_allow_html=True)

# Header
st.title("BioTech OS Pro")

# Sidebar
with st.sidebar:
    st.markdown("### SETTINGS")
    
    uploaded_file = st.file_uploader("Upload file", type=["fasta", "txt"])
    mol_type = st.selectbox("Sequence type", ["DNA", "RNA", "Protein"])
    
    st.markdown("---")
    st.markdown("### PARAMETERS")
    
    min_orf = st.slider("Min ORF length (bp)", 30, 300, 100)
    gc_window = st.slider("GC-skew window", 50, 500, 100, 50)
    
    selected_enzymes = st.multiselect(
        "Restriction enzymes",
        list(RESTRICTION_ENZYMES.keys()),
        default=['EcoRI', 'HindIII', 'BamHI']
    )
    
    st.markdown("---")
    st.markdown("### COMPARISON")
    compare_mode = st.checkbox("Enable sequence comparison")

# Main content
input_seq = ""
second_seq = ""

if uploaded_file:
    content = uploaded_file.getvalue().decode("utf-8")
    if content.startswith(">"):
        lines = content.split('\n')
        input_seq = "".join([line.strip() for line in lines[1:] if line and not line.startswith(">")])
    else:
        input_seq = content.strip()
else:
    default_seq = """ATGGCCATTGATGCGCTAGCTAGCTAGCTACGATCGATCGATCGATCG
    ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"""
    
    input_seq = st.text_area("Enter main sequence:", 
                           value=default_seq,
                           height=150).upper().replace(" ", "").replace("\n", "").strip()

if compare_mode:
    second_seq = st.text_area("Enter second sequence for comparison:",
                            height=100,
                            placeholder="ATGC... or MAIVMRS...").upper().replace(" ", "").replace("\n", "").strip()

if input_seq:
    # Validation
    valid_chars = {
        'DNA': set('ATCGN'),
        'RNA': set('AUCGN'),
        'Protein': set('ACDEFGHIKLMNPQRSTVWY*_')
    }
    
    current_set = valid_chars.get(mol_type, set())
    invalid_chars = set(input_seq) - current_set
    
    if invalid_chars and mol_type != 'Protein':
        st.error(f"Invalid characters: {', '.join(invalid_chars)}")
        st.stop()
    
    # Data processing
    if mol_type == "DNA":
        dna_seq = input_seq
        rna_seq = dna_seq.replace('T', 'U')
        protein_seq = "".join([GENETIC_CODE.get(dna_seq[i:i+3], "?") 
                              for i in range(0, len(dna_seq)-2, 3)])
        rev_comp = get_reverse_complement(dna_seq)
    elif mol_type == "RNA":
        rna_seq = input_seq
        dna_seq = rna_seq.replace('U', 'T')
        protein_seq = "".join([GENETIC_CODE.get(dna_seq[i:i+3], "?") 
                              for i in range(0, len(dna_seq)-2, 3)])
        rev_comp = get_reverse_complement(dna_seq)
    else:
        protein_seq = input_seq
        dna_seq = ""
        rna_seq = ""
    
    # Calculate pI
    pI = calculate_isoelectric_point_simple(protein_seq) if protein_seq else 7.0
    
    # METRICS
    col1, col2, col3, col4 = st.columns(4)
    
    if mol_type in ["DNA", "RNA"]:
        gc_content = ((dna_seq.count('G') + dna_seq.count('C')) / len(dna_seq) * 100 
                     if len(dna_seq) > 0 else 0)
        
        metrics = [
            ("Length", f"{len(dna_seq):,} bp"),
            ("GC Content", f"{gc_content:.1f}%"),
            ("Isoelectric Point", f"{pI:.2f}"),
            ("Protein Mass", f"{calculate_protein_mass(protein_seq):,.0f} Da")
        ]
    else:
        metrics = [
            ("Length", f"{len(protein_seq)} aa"),
            ("Isoelectric Point", f"{pI:.2f}"),
            ("Molecular Weight", f"{calculate_protein_mass(protein_seq):,.0f} Da"),
            ("Hydropathy", f"{np.mean([HYDROPATHY_INDEX.get(aa, 0) for aa in protein_seq]):.2f}")
        ]
    
    for i, (label, val) in enumerate(metrics):
        with [col1, col2, col3, col4][i]:
            st.markdown(f'''
                <div class="metric-box">
                    <small style="color: #8B949E;">{label}</small><br>
                    <span class="metric-val">{val}</span>
                </div>
            ''', unsafe_allow_html=True)
    
    # TABS
    tabs = st.tabs(["Sequences", "ORF Finder", "Restriction", 
                   "Isoelectric Point", "Plasmid Map", 
                   "Alignment", "Statistics", "Export"])
    
    # Tab 1: Sequences
    with tabs[0]:
        st.subheader("Sequences")
        
        col_seq1, col_seq2 = st.columns(2)
        
        with col_seq1:
            if mol_type in ["DNA", "RNA"]:
                st.write("**DNA (5'→3'):**")
                st.code(dna_seq[:500] + ("..." if len(dna_seq) > 500 else ""))
                
                st.write("**RNA transcript:**")
                st.code(rna_seq[:500] + ("..." if len(rna_seq) > 500 else ""))
        
        with col_seq2:
            if mol_type in ["DNA", "RNA"]:
                st.write("**Reverse complement:**")
                st.code(rev_comp[:500] + ("..." if len(rev_comp) > 500 else ""))
            
            st.write("**Protein sequence:**")
            st.code(protein_seq[:500] + ("..." if len(protein_seq) > 500 else ""))
    
    # Tab 2: ORF Finder
    with tabs[1]:
        if mol_type in ["DNA", "RNA"]:
            st.subheader("Open Reading Frames")
            
            orfs = find_orfs(dna_seq, min_orf)
            
            if orfs:
                st.write(f"**Found {len(orfs)} ORFs**")
                df_orfs = pd.DataFrame(orfs)
                st.dataframe(df_orfs[['Strand', 'Frame', 'Start', 'End', 'Length']], 
                           use_container_width=True)
                
                if len(orfs) > 0:
                    selected_orf = st.selectbox("ORF details:", 
                                              range(len(orfs)),
                                              format_func=lambda x: f"ORF {x+1}")
                    
                    with st.expander("Show details"):
                        st.write(f"**Protein:**")
                        st.code(orfs[selected_orf]['Protein'])
                        orf_pI = calculate_isoelectric_point_simple(orfs[selected_orf]['Protein'])
                        st.write(f"**pI:** {orf_pI:.2f}")
                        st.write(f"**Length:** {orfs[selected_orf]['Length']} bp")
            else:
                st.info("No ORFs found.")
        else:
            st.info("ORF analysis only for DNA/RNA.")
    
    # Tab 3: Restriction Analysis
    with tabs[2]:
        if mol_type == "DNA":
            st.subheader("Restriction Analysis")
            
            if selected_enzymes:
                enzymes = {k: RESTRICTION_ENZYMES[k] for k in selected_enzymes}
                sites = find_restriction_sites(dna_seq, enzymes)
                
                if sites:
                    col_map1, col_map2 = st.columns([2, 1])
                    
                    with col_map1:
                        st.write("**Restriction sites:**")
                        for site in sites:
                            st.write(f"• **{site['Enzyme']}**: {site['Position']} bp (sequence: {site['Sequence']})")
                    
                    with col_map2:
                        df_sites = pd.DataFrame(sites)
                        st.dataframe(df_sites[['Enzyme', 'Position', 'Sequence']], 
                                   use_container_width=True, height=300)
                        
                        # Fragments
                        if len(sites) >= 2:
                            positions = sorted([s['Position'] for s in sites])
                            positions.append(len(dna_seq))
                            positions.insert(0, 0)
                            
                            fragments = []
                            prev = 0
                            for pos in positions[1:]:
                                frag = pos - prev
                                if frag > 0:
                                    fragments.append(frag)
                                prev = pos
                            
                            st.write(f"**Fragments:** {len(fragments)}")
                            st.write(f"**Sizes (bp):** {', '.join(map(str, fragments))}")
                else:
                    st.info("No restriction sites found.")
            else:
                st.warning("Select enzymes for analysis.")
        else:
            st.info("Restriction analysis only for DNA.")
    
    # Tab 4: Isoelectric Point
    with tabs[3]:
        st.subheader("Isoelectric Point Analysis")
        
        if protein_seq:
            col_pi1, col_pi2 = st.columns(2)
            
            with col_pi1:
                # pI information
                acidic = protein_seq.count('D') + protein_seq.count('E')
                basic = protein_seq.count('R') + protein_seq.count('K') + protein_seq.count('H')
                
                st.markdown(f'''
                <div class="info-box">
                    <h4 style="color: #10B981; margin-top: 0;">Isoelectric Point = {pI:.2f}</h4>
                    <p style="color: #8B949E;">
                        pH where protein net charge is zero
                    </p>
                </div>
                ''', unsafe_allow_html=True)
                
                st.write("**Charged amino acids:**")
                st.metric("Acidic (D, E)", f"{acidic}")
                st.metric("Basic (R, K, H)", f"{basic}")
                
                # Interpretation
                if pI < 6:
                    st.success("**Acidic protein** - negative charge at physiological pH 7.4")
                elif pI > 8:
                    st.success("**Basic protein** - positive charge at physiological pH 7.4")
                else:
                    st.success("**Neutral protein** - charge near zero at pH 7.4")
            
            with col_pi2:
                # Charge vs pH plot
                st.markdown("<div class='blast-container'>", unsafe_allow_html=True)
                st.write("**Charge vs pH**")
                
                # Create simple text plot
                plot_text = create_simple_charge_plot(protein_seq)
                st.markdown(f'<div class="simple-plot">{plot_text}</div>', unsafe_allow_html=True)
                
                st.markdown("</div>", unsafe_allow_html=True)
        else:
            st.info("No protein sequence for pI analysis.")
    
    # Tab 5: Plasmid Map
    with tabs[4]:
        if mol_type == "DNA":
            st.subheader("Plasmid Map")
            
            enzymes = {k: RESTRICTION_ENZYMES[k] for k in selected_enzymes}
            sites = find_restriction_sites(dna_seq, enzymes)
            
            if sites:
                # Text plasmid map
                map_text = create_text_plasmid_map(dna_seq, sites)
                st.markdown(f'<div class="simple-plot">{map_text}</div>', unsafe_allow_html=True)
                
                # Sites table
                st.write("**Restriction sites:**")
                df_sites = pd.DataFrame(sites)
                st.dataframe(df_sites[['Enzyme', 'Position', 'Sequence', 'Cut Site']], 
                           use_container_width=True)
            else:
                st.info("No restriction sites found for plasmid map.")
        else:
            st.info("Plasmid map only for DNA.")
    
    # Tab 6: Alignment
    with tabs[5]:
        st.subheader("Sequence Alignment")
        
        if compare_mode and second_seq:
            # Check second sequence
            invalid_chars_2 = set(second_seq) - current_set
            if invalid_chars_2 and mol_type != 'Protein':
                st.error(f"Invalid characters in second sequence: {', '.join(invalid_chars_2)}")
            else:
                # Perform alignment
                alignment = needleman_wunsch(input_seq, second_seq)
                
                col_align1, col_align2 = st.columns(2)
                
                with col_align1:
                    st.write("**Alignment results:**")
                    st.metric("Score", f"{alignment['score']:.1f}")
                    st.metric("Identity", f"{alignment['identity']:.1f}%")
                    st.metric("Matches", alignment['matches'])
                    st.metric("Gaps", alignment['gaps'])
                
                with col_align2:
                    st.write("**Information:**")
                    st.write(f"• Sequence 1: {len(input_seq)} characters")
                    st.write(f"• Sequence 2: {len(second_seq)} characters")
                    st.write(f"• Alignment length: {len(alignment['alignment1'])}")
                
                # Display alignment
                st.markdown("---")
                st.write("**Aligned sequences (first 200 characters):**")
                
                # Show first 200 characters
                align1 = alignment['alignment1'][:200]
                align2 = alignment['alignment2'][:200]
                
                html_output = ""
                for i in range(0, len(align1), 50):
                    block1 = align1[i:i+50]
                    block2 = align2[i:i+50]
                    
                    # Highlighting
                    line1 = ""
                    line2 = ""
                    matches = ""
                    
                    for a, b in zip(block1, block2):
                        if a == b and a != '-':
                            line1 += f'<span class="alignment-match">{a}</span>'
                            line2 += f'<span class="alignment-match">{b}</span>'
                            matches += "|"
                        elif a == '-' or b == '-':
                            line1 += f'<span class="alignment-gap">{a}</span>'
                            line2 += f'<span class="alignment-gap">{b}</span>'
                            matches += " "
                        else:
                            line1 += f'<span class="alignment-mismatch">{a}</span>'
                            line2 += f'<span class="alignment-mismatch">{b}</span>'
                            matches += "."
                    
                    html_output += f"""
                    <div style="font-family: 'Courier New', monospace; background: #1E2229; 
                                padding: 10px; margin: 5px 0; border-radius: 5px;">
                        <div>Seq1: {line1}</div>
                        <div style="color: #10B981; padding-left: 6px;">{matches}</div>
                        <div>Seq2: {line2}</div>
                    </div>
                    """
                
                st.markdown(html_output, unsafe_allow_html=True)
        else:
            if not compare_mode:
                st.info("Enable 'Sequence comparison' in sidebar.")
            else:
                st.info("Enter second sequence for comparison.")
    
    # Tab 7: Statistics
    with tabs[6]:
        st.subheader("Statistical Analysis")
        
        if mol_type in ["DNA", "RNA"]:
            col_stat1, col_stat2 = st.columns(2)
            
            with col_stat1:
                # Nucleotide composition
                st.write("**Nucleotide composition:**")
                bases = ['A', 'T', 'C', 'G'] if mol_type == 'DNA' else ['A', 'U', 'C', 'G']
                counts = {base: dna_seq.count(base) for base in bases}
                
                df_bases = pd.DataFrame({
                    'Nucleotide': bases,
                    'Count': [counts[b] for b in bases],
                    'Percent': [f"{counts[b]/len(dna_seq)*100:.1f}%" for b in bases]
                })
                st.dataframe(df_bases, use_container_width=True)
            
            with col_stat2:
                # GC-skew
                st.write("**GC-skew analysis:**")
                positions, skews = calculate_gc_skew(dna_seq, gc_window)
                
                if positions:
                    # Simple line chart
                    skew_df = pd.DataFrame({'Position': positions, 'GC-skew': skews})
                    st.line_chart(skew_df.set_index('Position'))
                
                # Hydropathy
                if protein_seq and len(protein_seq) > 10:
                    st.write("**Hydropathy profile:**")
                    positions, values = calculate_hydropathy_profile(protein_seq)
                    
                    if positions:
                        hyd_df = pd.DataFrame({'Position': positions, 'Hydropathy': values})
                        st.line_chart(hyd_df.set_index('Position'))
        
        # Amino acid composition
        if protein_seq:
            st.write("**Amino acid composition:**")
            aa_counts = Counter(protein_seq)
            aa_data = []
            
            for aa, count in aa_counts.items():
                freq = count / len(protein_seq) * 100
                weight = AA_WEIGHTS.get(aa, 0)
                hydropathy = HYDROPATHY_INDEX.get(aa, 'N/A')
                aa_data.append({
                    'AA': aa,
                    'Count': count,
                    'Frequency (%)': f"{freq:.1f}%",
                    'Weight': f"{weight:.1f}",
                    'Hydropathy': hydropathy
                })
            
            df_aa = pd.DataFrame(aa_data)
            st.dataframe(df_aa, use_container_width=True)
    
    # Tab 8: EXPORT
    with tabs[7]:
        st.subheader("Data Export")
        
        st.markdown("<div class='export-section'>", unsafe_allow_html=True)
        
        col_exp1, col_exp2 = st.columns(2)
        
        with col_exp1:
            st.write("**Export Sequences**")
            
            # FASTA export
            if mol_type in ["DNA", "RNA"]:
                fasta_content = f">Original_DNA_sequence\n{dna_seq}\n\n>RNA_transcript\n{rna_seq}\n\n>Protein_translation\n{protein_seq}\n\n>Reverse_complement\n{rev_comp}"
                
                st.download_button(
                    label="Download all sequences (FASTA)",
                    data=fasta_content,
                    file_name="all_sequences.fasta",
                    mime="text/plain",
                    use_container_width=True
                )
                
                st.download_button(
                    label="Download DNA only (FASTA)",
                    data=f">DNA_sequence\n{dna_seq}",
                    file_name="dna_sequence.fasta",
                    mime="text/plain",
                    use_container_width=True
                )
            
            if protein_seq:
                st.download_button(
                    label="Download protein sequence (FASTA)",
                    data=f">Protein_sequence\n{protein_seq}",
                    file_name="protein_sequence.fasta",
                    mime="text/plain",
                    use_container_width=True
                )
        
        with col_exp2:
            st.write("**Export Analysis Results**")
            
            # Prepare data for CSV
            analysis_data = []
            
            if mol_type in ["DNA", "RNA"]:
                analysis_data.append({"Parameter": "Type", "Value": mol_type})
                analysis_data.append({"Parameter": "Length", "Value": f"{len(dna_seq)} bp"})
                analysis_data.append({"Parameter": "GC Content", "Value": f"{gc_content:.1f}%"})
                analysis_data.append({"Parameter": "Isoelectric Point", "Value": f"{pI:.2f}"})
                analysis_data.append({"Parameter": "Protein Mass", "Value": f"{calculate_protein_mass(protein_seq):,.0f} Da"})
                
                # Nucleotide composition
                bases = ['A', 'T', 'C', 'G'] if mol_type == 'DNA' else ['A', 'U', 'C', 'G']
                for base in bases:
                    count = dna_seq.count(base)
                    freq = count / len(dna_seq) * 100
                    analysis_data.append({"Parameter": f"Nucleotide {base}", "Value": f"{count} ({freq:.1f}%)"})
            
            # Convert to CSV
            df_analysis = pd.DataFrame(analysis_data)
            csv_data = df_analysis.to_csv(index=False)
            
            st.download_button(
                label="Download analysis (CSV)",
                data=csv_data,
                file_name="analysis_results.csv",
                mime="text/csv",
                use_container_width=True
            )
            
            # JSON export
            json_results = {
                "sequence_type": mol_type,
                "length": len(dna_seq) if mol_type in ["DNA", "RNA"] else len(protein_seq),
                "gc_content": f"{gc_content:.1f}%" if mol_type in ["DNA", "RNA"] else None,
                "isoelectric_point": f"{pI:.2f}",
                "protein_molecular_weight": f"{calculate_protein_mass(protein_seq):,.0f} Da",
                "protein_sequence": protein_seq if protein_seq else None
            }
            
            if mol_type in ["DNA", "RNA"]:
                json_results["dna_sequence"] = dna_seq
                json_results["rna_sequence"] = rna_seq
                json_results["reverse_complement"] = rev_comp
            
            json_data = json.dumps(json_results, indent=2)
            
            st.download_button(
                label="Download full report (JSON)",
                data=json_data,
                file_name="full_analysis.json",
                mime="application/json",
                use_container_width=True
            )
        
        st.markdown("</div>", unsafe_allow_html=True)
        
        # NCBI BLAST section
        st.markdown("---")
        st.write("**NCBI BLAST Analysis**")
        
        col_blast1, col_blast2 = st.columns(2)
        
        with col_blast1:
            # Form for BLAST
            with st.form("blast_form"):
                st.write("Configure BLAST parameters:")
                
                if mol_type in ["DNA", "RNA"]:
                    blast_program = st.selectbox(
                        "BLAST program",
                        ["blastn", "blastp", "blastx", "tblastn", "tblastx"],
                        help="Select appropriate BLAST program for your sequence"
                    )
                else:
                    blast_program = "blastp"
                
                database = st.selectbox(
                    "Database",
                    ["nr", "refseq_rna", "refseq_genomic", "pdb"] if mol_type in ["DNA", "RNA"] else ["nr", "pdb", "swissprot"],
                    help="Target database for search"
                )
                
                e_value = st.number_input("E-value threshold", min_value=1e-100, max_value=10.0, value=0.001, format="%.3f")
                
                submitted = st.form_submit_button("Open NCBI BLAST")
                
                if submitted:
                    # Encode sequence for URL
                    sequence = dna_seq if mol_type in ["DNA", "RNA"] else protein_seq
                    encoded_seq = urllib.parse.quote(sequence[:3000])  # Limit length
                    
                    # Form BLAST URL
                    blast_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM={blast_program}&PAGE_TYPE=BlastSearch&BLAST_SPEC=&QUERY={encoded_seq}&DATABASE={database}&EXPECT={e_value}"
                    
                    st.success("BLAST query created!")
                    st.markdown(f"""
                    <div class="blast-container">
                        <h4>BLAST analysis link:</h4>
                        <p><a href="{blast_url}" target="_blank">{blast_url[:100]}...</a></p>
                        <p><strong>Parameters:</strong></p>
                        <ul>
                            <li>Program: {blast_program}</li>
                            <li>Database: {database}</li>
                            <li>E-value: {e_value}</li>
                            <li>Sequence length: {len(sequence)} characters</li>
                        </ul>
                    </div>
                    """, unsafe_allow_html=True)
        
        with col_blast2:
            st.write("**BLAST Information:**")
            st.info("""
            **What is BLAST?**
            
            BLAST (Basic Local Alignment Search Tool) is an algorithm for comparing 
            biological sequences such as DNA, RNA or proteins.
            
            **Recommendations:**
            - For DNA use **blastn**
            - For proteins use **blastp**
            - For short sequences increase E-value
            - **nr** - most complete database
            - **pdb** - 3D structure database
            """)
            
            # Quick links
            st.write("**Quick links:**")
            
            if mol_type in ["DNA", "RNA"]:
                dna_blast_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&QUERY={urllib.parse.quote(dna_seq[:3000])}"
                st.markdown(f'<a href="{dna_blast_url}" target="_blank"><button style="width: 100%; margin: 5px 0;">BLAST for DNA</button></a>', unsafe_allow_html=True)
            
            if protein_seq:
                protein_blast_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&QUERY={urllib.parse.quote(protein_seq[:3000])}"
                st.markdown(f'<a href="{protein_blast_url}" target="_blank"><button style="width: 100%; margin: 5px 0;">BLAST for protein</button></a>', unsafe_allow_html=True)

else:
    # Start screen
    st.markdown("""
    <div style='text-align: center; padding: 50px;'>
        <h2 style='color: #10B981;'>BioTech OS Pro</h2>
        <p style='color: #8B949E; font-size: 1.1rem; margin-bottom: 40px;'>
            Professional bioinformatics analyzer with full functionality
        </p>
        
        <div style='display: flex; justify-content: center; flex-wrap: wrap; gap: 20px;'>
            <div style='background-color: #1E2229; padding: 20px; border-radius: 10px; width: 200px;'>
                <h4 style='color: #10B981;'>Isoelectric Point</h4>
                <p style='color: #8B949E; font-size: 0.9rem;'>Calculate protein pI</p>
            </div>
            <div style='background-color: #1E2229; padding: 20px; border-radius: 10px; width: 200px;'>
                <h4 style='color: #10B981;'>Plasmid Map</h4>
                <p style='color: #8B949E; font-size: 0.9rem;'>Restriction site visualization</p>
            </div>
            <div style='background-color: #1E2229; padding: 20px; border-radius: 10px; width: 200px;'>
                <h4 style='color: #10B981;'>Data Export</h4>
                <p style='color: #8B949E; font-size: 0.9rem;'>FASTA, CSV, JSON formats</p>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)

# Footer
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #8B949E; font-size: 0.9rem;'>
    BioTech OS Pro v3.1 | Full functionality with working export and BLAST
</div>
""", unsafe_allow_html=True)
