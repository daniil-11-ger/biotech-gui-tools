import tkinter as tk
from tkinter import messagebox

def calculate_tm():
    """Рассчитывает температуру отжига праймера (Tm)"""
    seq = seq_entry.get().upper().strip()
    
    # Проверка на корректность ДНК
    if not seq or any(base not in 'ATGC' for base in seq):
        messagebox.showerror("Ошибка", "Введите корректную последовательность (A, T, G, C)")
        return

    # Формула Уоллеса (для коротких праймеров)
    # Tm = 2*(A+T) + 4*(G+C)
    a_t = seq.count('A') + seq.count('T')
    g_c = seq.count('G') + seq.count('C')
    tm = 2 * a_t + 4 * g_c
    
    gc_content = (g_c / len(seq)) * 100
    
    result_text = f"Последовательность: {seq}\n" \
                  f"Длина: {len(seq)} bp\n" \
                  f"GC-состав: {gc_content:.1f}%\n" \
                  f"Tm (примерная): {tm}°C"
    
    messagebox.showinfo("Результат анализа", result_text)

# Настройка главного окна
window = tk.Tk()
window.title("Biotech Lab Assistant")
window.geometry("350x300")
window.configure(bg='#F0F8FF') # Светло-голубой "лабораторный" фон

# Оформление
header = tk.Label(window, text="Primer Analysis Tool", font=("Helvetica", 16, "bold"), 
                  bg='#F0F8FF', fg='#2F4F4F')
header.pack(pady=15)

instruction = tk.Label(window, text="Введите последовательность праймера (5'->3'):", 
                       bg='#F0F8FF', font=("Helvetica", 10))
instruction.pack()

seq_entry = tk.Entry(window, font=("Courier", 12), width=25, justify='center')
seq_entry.pack(pady=10)

# Кнопка расчета
calc_btn = tk.Button(window, text="Рассчитать Tm и GC", 
                     command=calculate_tm, 
                     bg='#4682B4', fg='white', 
                     font=("Helvetica", 10, "bold"), 
                     padx=10, pady=5)
calc_btn.pack(pady=20)

# Инфо-блок внизу
footer = tk.Label(window, text="BioTech 3rd Year Project", 
                  bg='#F0F8FF', fg='#778899', font=("Helvetica", 8, "italic"))
footer.pack(side="bottom", pady=10)

window.mainloop()
