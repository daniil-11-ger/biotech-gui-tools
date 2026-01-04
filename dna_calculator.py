import streamlit as st

# --- Bioinformatics Logic ---

def calculate_dna(seq):
    at_count = seq.count('A') + seq.count('T')
    gc_count = seq.count('G') + seq.count('C')
    tm = 2 * at_count + 4 * gc_count
    gc_percent = (gc_count / len(seq)) * 100 if len(seq) > 0 else 0
    return tm, gc_percent

def transcribe_dna(seq):
    return seq.replace('T', 'U')

# --- Streamlit Web Interface ---

st.set_page_config(page_title="BioTech Lab Assistant v1.1", page_icon="🧬")

st.title("🧬 DNA Analysis Tool")
st.markdown("Developed for Biotechnology Research | 2026")

# Input Field
seq = st.text_input("Enter DNA Sequence (5'->3'):", placeholder="Example: ATGCGTAC...").upper().strip()

# Buttons (Action Columns)
col1, col2 = st.columns(2)

if seq:
    # Validation
    if any(base not in 'ATGC' for base in seq):
        st.error("Please enter a valid DNA sequence (A, T, G, C)")
    else:
        with col1:
            if st.button("Calculate Tm & GC%"):
                tm, gc_percent = calculate_dna(seq)
                st.success(f"**Results:**\n\nLength: {len(seq)} bp\n\nGC-Content: {gc_percent:.1f}%\n\nMelting Temp (Tm): {tm}°C")
        
        with col2:
            if st.button("Transcription (DNA->RNA)"):
                rna = transcribe_dna(seq)
                st.info(f"**mRNA Sequence:**\n\n{rna}")
else:
    st.write("Enter a sequence to see analysis options.")
