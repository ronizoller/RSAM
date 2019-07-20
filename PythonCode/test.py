import tkinter as tk
from decimal import *
import RSAMfinder
from tkinter import ttk
import threading

patterns = ['p1', 'p2']

class Main_Frame(object):
    def __init__(self, top=None):
        # save root reference
        self.top = top
        # set title bar
        self.top.title("RSAM-finder")

        # start button calls the "initialization" function bar_init, you can pass a variable in here if desired
        labels = [['Specie tree extension', 'string', None], ['k', 'int', '50'],
                  ['Threshold Edges\nin Subtree', 'float', '0.1'],
                  ['HT event cost', 'float', '1'], ['Duplication event Cost', 'float', '1'],
                  ['Speciation event Cost', 'float', '0'],['Loss event Cost', 'float', '0'], ['Gamma', 'float', '1'],
                  ['p', 'float', '0.05'], ['Number of Vertices to find', 'int', 1]]

        entries = makeform(top, labels)
        v = tk.BooleanVar()
        MODES = [("Double-Mode", True), ("Single-Mod", False)]
        for text, mode in MODES:
            b = tk.Radiobutton(top, indicatoron=0, text=text, variable=v, value=mode)
            b.pack(side=tk.LEFT, padx=5, pady=5)

        self.b1 = tk.Button(top, text='Run',
                            command=lambda: self.bar_init(2500,entries,top,v))
        self.b1.pack(side=tk.LEFT, padx=5, pady=5)
        b2 = tk.Button(top, text='Quit', command=top.quit)
        b2.pack(side=tk.LEFT, padx=5, pady=5)

        # run mainloop
        self.top.mainloop()

    def bar_init(self, var, ent,top,v):
        self.window = tk.Toplevel(root)
        top.title("RSAM-finder")

        msg = tk.Message(self.window, text='RSAM-finder in progress...')
        msg.grid(row=0,column=0)
        self.load_bar = ttk.Progressbar(self.window)
        self.load_bar.grid(row=1,column=0)
        # first layer of isolation, note var being passed along to the self.start_bar function
        # target is the function being started on a new thread, so the "bar handler" thread
        self.start_bar_thread = threading.Thread(target=self.start_bar, args=(var,ent,v,))
        # start the bar handling thread
        self.start_bar_thread.start()

    def start_bar(self, var, ent,v):
        # the load_bar needs to be configured for indeterminate amount of bouncing
        self.load_bar.config(mode='determinate', maximum=100, value=0)
        # 8 here is for speed of bounce
        self.load_bar.start(8)
        # start the work-intensive thread, again a var can be passed in here too if desired
        self.work_thread = threading.Thread(target=self.work_task, args=(var,ent,v,))
        self.work_thread.start()
        # close the work thread
        self.work_thread.join()
        # stop the indeterminate bouncing
        self.load_bar.stop()
        # reconfigure the bar so it appears reset
        self.load_bar.config(value=0, maximum=0)
        self.window.destroy()
        msg = tk.Label(self.top, width=30,text='Results can be found in /data/results.')
        msg.pack()
        self.b1.destroy()


    def work_task(self,var,entries,double):
        p1_EV = []
        p1_colors = None
        p1_dis = False
        p2_EV = []
        p2_colors = None
        p2_dis = False
        for e in entries:
            if e[1] != None:
                for p in patterns:
                    if e[0].find(p) != -1:
                        if e[0][3:].find('event') != -1:
                            eval(p + '_EV').append(e[1].get())
                        if e[0][3:].find('colors') != -1:
                            exec(p + '_colors = e[1].get()')
                        if e[0][3:].find('dis') != -1:
                            exec(p + '_dis = e[1].get()')
        p1 = (p1_EV, p1_colors, p1_dis)
        if not double:
            p2 = (p2_EV, p2_colors, p2_dis)
        else:
            p2 = (None,None,False)
        vars = []
        for entry in entries:
            if entry[0].find('p1') == -1 and entry[0].find('p2') == -1:
                text = entry[1].get()
                vars.append(text)
        speciesTreespecification, k, TH_edges, HT_cost, D_cost, S_cost,loss_cost, gamma, p, number_of_planted_vertices,create_sigma = vars
        RSAMfinder.main(speciesTreespecification, int(k), Decimal(TH_edges), int(HT_cost), int(D_cost), int(S_cost),int(loss_cost),
                        Decimal(gamma), Decimal(p), int(number_of_planted_vertices), p1, p2, True,create_sigma)

def makeform(root, fields):
    def intValidation(S):
        flag = True
        for s in S:
            if s not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                flag = False
        if flag:
            return True
        t3.bell()  # .bell() plays that ding sound telling you there was invalid input
        return False
    def floatValidation(S):
        flag = True
        for s in S:
            if s not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9','.']:
                flag = False
        if flag:
            return True
        t3.bell()  # .bell() plays that ding sound telling you there was invalid input
        return False

    entries = []
    t3 = tk.Toplevel(root)
    intVal = (t3.register(intValidation), '%S')
    floatVal = (t3.register(floatValidation), '%S')
    lab = tk.Label(root, text='RSAM-finder')
    lab.config(width='50')
    lab.pack()
    for field in fields:
        row = tk.Frame(root)
        lab = tk.Label(row, width=30, text=field[0])
        if field[1] == 'float':
            ent = tk.Entry(row, validate='key',vcmd=floatVal)
        elif field[1] == 'int':
            ent = tk.Entry(row, validate='key', vcmd=intVal)
        else:
            ent = tk.Entry(row, validate='key')
        if field[2] != None:
            ent.insert("end",field[2])
        row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
        lab.pack(side=tk.LEFT)
        ent.pack(side=tk.RIGHT, expand=tk.YES)
        entries.append((field[0], ent))
    for p in patterns:
        row = tk.Frame(root)
        names = ['Event1','Event2','Event3','Colors','Distance']
        for name in names:
            lab = tk.Label(row, text=name,width=10)
            lab.grid(row=0, column=names.index(name)+1)
        lab = tk.Label(row, width=10, text=p)
        lab.grid(row=1,column=0)
        number_of_events = 3
        for i in range(0,number_of_events):
            v = tk.StringVar(root)
            EV = ['S','D','HT']
            EV = tk.OptionMenu(row, v, *EV)
            EV.grid(row=1, column=i+1)
            entries.append((p+'_event'+str(i),v))
        colors = ['red','black','None']
        v = tk.StringVar(root)
        entries.append((p + 'colors', v))
        colors_op = tk.OptionMenu(row, v, *colors)
        colors_op.config(width=10)
        colors_op.grid(row=1, column=number_of_events+1)
        v = tk.BooleanVar(root)
        c = tk.Checkbutton(row,variable=v)
        entries.append((p + '_dis', v))
        c.grid(row=1, column=number_of_events+2)
        row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
    row = tk.Frame(root)
    lab = tk.Label(row, text='create sigma from FASTA file', width=40)
    lab.grid(row=0, column=0)
    create_sigma = tk.BooleanVar(root)
    c = tk.Checkbutton(row, variable=create_sigma)
    c.grid(row=0, column=1)
    entries.append(('create_sigma', create_sigma))
    row.pack(side=tk.TOP, fill=tk.X, padx=0, pady=5)
    return entries

if __name__ == '__main__':
    root = tk.Tk()
    Main_Frame(root)
