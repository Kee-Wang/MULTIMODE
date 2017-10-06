import tkinter as tk

root = tk.Tk()
mylist = ['a','b','c','d','e']

for i, x in enumerate(mylist):
    label = tk.Label(root, text="Label "+str(i))
    label.grid(row=i+1, column=1)
    label.bind("<Enter>", lambda e, x=x: e.widget.config(text=x))
    label.bind("<Leave>", lambda e, i=i: e.widget.config(text="Label "+str(i)))

root.mainloop()