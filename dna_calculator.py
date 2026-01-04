import streamlit as st

# --- Логика расчетов (твоя база) ---
def calculate_dna(seq):
    at_count = seq.count('A') + seq.count('T')
    gc_count = seq.count('G') + seq.count('C')
    tm = 2 * at_count + 4 * gc_count
    gc_percent = (gc_count / len(seq)) * 100 if len(seq) > 0 else 0
    return tm, gc_percent

def transcribe_dna(seq):
    return seq.replace('T', 'U')

# --- Настройка веб-интерфейса Streamlit ---
st.set_page_config(page_title="BioTech Lab Assistant", page_icon="🧬")

st.title("🧬 BioTech Lab Assistant")
st.markdown("Professional Laboratory DNA Analysis Tool")

# Поле ввода
seq_input = st.text_input("Enter DNA Sequence (5'->3'):", placeholder="ATGC...").upper().strip()

if seq_input:
    # Проверка на ошибки
    if any(base not in 'ATGC' for base in seq_input):
        st.error("Invalid sequence! Please use only A, T, G, C.")
    else:
        # Кнопки действий
        col1, col2 = st.columns(2)
        
        with col1:
            if st.button("Calculate Tm & GC%"):
                tm, gc_val = calculate_dna(seq_input)
                st.success(f"**Tm:** {tm}°C  \n**GC:** {gc_val:.1f}%")
        
        with col2:
            if st.button("Transcription"):
                rna = transcribe_dna(seq_input)
                st.info(f"**mRNA:** \n{rna}")
        
        # Маленький бонус: визуализация состава
        st.bar_chart({"A": seq_input.count('A'), "T": seq_input.count('T'), 
                      "G": seq_input.count('G'), "C": seq_input.count('C')})
else:
    st.info("Waiting for sequence input...")
