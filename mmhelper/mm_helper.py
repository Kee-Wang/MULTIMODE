from tkinter import *


class AllWidgets():
    def __init__(self, root):
        self.frame = Frame(root)

        self.build_window()

        self.frame.pack()

        menubar = Menu(root)
        root['menu'] = menubar

        menu_file = Menu(menubar)
        menu_file.add_command(label="Quit", command=self.quit)

        menubar.add_cascade(menu=menu_file, label='File')

    def build_window(self):

        Label(self.frame, text="MULTIMODE's little helper").pack(side=TOP)
        Button(self.frame, text="Check VCI size now").pack(side=TOP)

        Entry(self.frame, text="Entry").pack(side=TOP)
        lb = Listbox(self.frame, height=3)
        lb.pack(side=TOP)
        for x in ('Listbox 1', 'Listbox 2', 'Listbox 3'):
            lb.insert('end', x)
        # RadioButton(self.frame, text="Radio Button").pack(side=TOP)
        LabelFrame(self.frame, text="Label Frame").pack(side=TOP)
        Scale(self.frame, from_=0, to=200, orient=HORIZONTAL).pack(side=TOP)
        Scrollbar(self.frame, orient=HORIZONTAL).pack(side=TOP, fill=X)
        Spinbox(self.frame, from_=0, to=200).pack(side=TOP)
        p = PhotoImage(file="1.gif")
        l = Label(self.frame, image=p)
        l.image = p
        l.pack(side=TOP)

        t = Text(self.frame, height=5, width=20)
        t.pack(side=TOP)
        t.insert(END, "User-defined data here. \n \nThis part will be read line-by-line in `user.XXX.f`")

    def quit(self):
        self.frame.quit()


if __name__ == '__main__':
    root = Tk()

#    app = AllWidgets(root)

 #   root.mainloop()
root = Tk() #root means to creat a new window
theLabel = Label(root,text="MSA PES fitting package").pack(side=TOP)
Entry(text="Entry").pack(side=TOP)
#theLabel.pack() #label pops out really fast and close really fast.

root.mainloop() #This means have window continuous on the screen until
