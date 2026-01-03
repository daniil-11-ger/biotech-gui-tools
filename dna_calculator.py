import tkinter as tk
from tkinter import messagebox

# --- Bioinformatics Logic ---

def calculate_dna():
    """Calculates Tm and GC-content of the sequence"""
    seq = seq_entry.get().upper().strip()
    
    # Validation: Ensure only nucleotide bases are entered
    if not seq or any(base not in 'ATGC' for base in seq):
        messagebox.showerror("Input Error", "Please enter a valid DNA sequence (A, T, G, C)")
        return

    # Wallace Rule Calculation: Tm = 2*(A+T) + 4*(G+C)
    at_count = seq.count('A') + seq.count('T')
    gc_count = seq.count('G') + seq.count('C')
    tm = 2 * at_count + 4 * gc_count
    
    # GC-Content Calculation
    gc_percent = (gc_count / len(seq)) * 100
    
    result = f"Analysis for: {seq}\n" \
             f"--------------------------\n" \
             f"Length: {len(seq)} bp\n" \
             f"GC-Content: {gc_percent:.1f}%\n" \
             f"Melting Temp (Tm): {tm}°C"
    
    messagebox.showinfo("Analysis Results", result)

def transcribe_dna():
    """Performs DNA to RNA transcription"""
    seq = seq_entry.get().upper().strip()
    if not seq or any(base not in 'ATGC' for base in seq):
        messagebox.showerror("Error", "Invalid sequence for transcription")
        return
    
    rna = seq.replace('T', 'U')
    messagebox.showinfo("Transcription", f"mRNA Sequence:\n\n{rna}")

# --- GUI Configuration ---

root = tk.Tk()
root.title("BioTech Lab Assistant v1.1")
root.geometry("400x350")
root.resizable(False, False)
root.configure(bg='#F5F5F5') # Professional lab-style background

# Header
header = tk.Label(root, text="DNA Analysis Tool", font=("Segoe UI", 16, "bold"), 
                  bg='#F5F5F5', fg='#2C3E50')
header.pack(pady=20)

# Input Field
input_label = tk.Label(root, text="Enter DNA Sequence (5'->3'):", 
                       bg='#F5F5F5', font=("Segoe UI", 10))
input_label.pack()

seq_entry = tk.Entry(root, font=("Consolas", 12), width=30, justify='center', 
                     fg='#2C3E50', relief="flat", highlightthickness=1)
seq_entry.pack(pady=10)

# Hint
hint = tk.Label(root, text="Example: ATGCGTAC... ", bg='#F5F5F5', 
                fg='#7F8C8D', font=("Segoe UI", 8))
hint.pack()

# Action Buttons
btn_frame = tk.Frame(root, bg='#F5F5F5')
btn_frame.pack(pady=20)

# Button 1: Tm/GC Analysis
calc_btn = tk.Button(btn_frame, text="Calculate Tm & GC%", command=calculate_dna, 
                     bg='#3498DB', fg='white', font=("Segoe UI", 10, "bold"), 
                     width=25, relief="flat", cursor="hand2")
calc_btn.pack(pady=5)

# Button 2: Transcription
trans_btn = tk.Button(btn_frame, text="Transcription (DNA->RNA)", command=transcribe_dna, 
                      bg='#2ECC71', fg='white', font=("Segoe UI", 10, "bold"), 
                      width=25, relief="flat", cursor="hand2")
trans_btn.pack(pady=5)

# Footer
footer = tk.Label(root, text="Developed for Biotechnology Research | 2026", 
                  bg='#F5F5F5', fg='#BDC3C7', font=("Segoe UI", 8))
footer.pack(side="bottom", pady=10)

root.mainloop()
