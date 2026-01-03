import tkinter as tk
from tkinter import messagebox

# --- Логика расчетов ---

def calculate_dna():
    """Рассчитывает Tm и GC-состав последовательности"""
    seq = seq_entry.get().upper().strip()
    
    # Валидация: проверка, что введены только нуклеотиды
    if not seq or any(base not in 'ATGC' for base in seq):
        messagebox.showerror("Ошибка ввода", "Введите корректную последовательность ДНК (A, T, G, C)")
        return

    # Расчет Tm по формуле Уоллеса: Tm = 2*(A+T) + 4*(G+C)
    at_count = seq.count('A') + seq.count('T')
    gc_count = seq.count('G') + seq.count('C')
    tm = 2 * at_count + 4 * gc_count
    
    # Расчет GC-состава
    gc_percent = (gc_count / len(seq)) * 100
    
    result = f"Анализ для: {seq}\n" \
             f"--------------------------\n" \
             f"Длина: {len(seq)} bp\n" \
             f"GC-состав: {gc_percent:.1f}%\n" \
             f"Температура отжига (Tm): {tm}°C"
    
    messagebox.showinfo("Результаты анализа", result)

def transcribe_dna():
    """Выполняет транскрипцию ДНК в РНК"""
    seq = seq_entry.get().upper().strip()
    if not seq or any(base not in 'ATGC' for base in seq):
        messagebox.showerror("Ошибка", "Некорректная последовательность для транскрипции")
        return
    
    rna = seq.replace('T', 'U')
    messagebox.showinfo("Транскрипция", f"Последовательность мРНК:\n\n{rna}")

# --- Настройка графического интерфейса ---

root = tk.Tk()
root.title("BioTech Lab Assistant v1.1")
root.geometry("400x350")
root.resizable(False, False)
root.configure(bg='#F5F5F5') # Нейтральный лабораторный фон

# Заголовок
header = tk.Label(root, text="DNA Analysis Tool", font=("Segoe UI", 16, "bold"), 
                  bg='#F5F5F5', fg='#2C3E50')
header.pack(pady=20)

# Поле ввода
input_label = tk.Label(root, text="Введите последовательность ДНК (5'->3'):", 
                       bg='#F5F5F5', font=("Segoe UI", 10))
input_label.pack()

seq_entry = tk.Entry(root, font=("Consolas", 12), width=30, justify='center', 
                     fg='#2C3E50', relief="flat", highlightthickness=1)
seq_entry.pack(pady=10)

# Подсказка
hint = tk.Label(root, text="Пример: ATGCGTAC... ", bg='#F5F5F5', 
                fg='#7F8C8D', font=("Segoe UI", 8))
hint.pack()

# Кнопки действий
btn_frame = tk.Frame(root, bg='#F5F5F5')
btn_frame.pack(pady=20)

# Кнопка 1: Анализ Tm/GC
calc_btn = tk.Button(btn_frame, text="Расчет Tm & GC%", command=calculate_dna, 
                     bg='#3498DB', fg='white', font=("Segoe UI", 10, "bold"), 
                     width=20, relief="flat", cursor="hand2")
calc_btn.pack(pady=5)

# Кнопка 2: Транскрипция
trans_btn = tk.Button(btn_frame, text="Транскрипция (DNA->RNA)", command=transcribe_dna, 
                      bg='#2ECC71', fg='white', font=("Segoe UI", 10, "bold"), 
                      width=20, relief="flat", cursor="hand2")
trans_btn.pack(pady=5)

# Подвал
footer = tk.Label(root, text="Designed for Biotechnology Students | 2026", 
                  bg='#F5F5F5', fg='#BDC3C7', font=("Segoe UI", 8))
footer.pack(side="bottom", pady=10)

root.mainloop()
