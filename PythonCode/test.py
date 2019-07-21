import tkinter as tk
from decimal import *
import RSAMfinder
from tkinter import ttk
import threading
from PIL import ImageTk, Image


class parameters_frame(object):
    def __init__(self,top,labels,specification):
        self.top = top
        self.v = None

        if specification == 'general':
            self.entries = makeform(self.top, labels)
        elif specification == 'patterns':
            self.entries = make_patterns_tab(self.top, ['p1','p2'])
        elif specification == 'figure_params':
            self.entries = make_figure_parameters(self,self.top)

    def work_task (self):
        return self.entries

    def callback_color(self,vars):
        print('TODO')

def make_patterns_tab(root, patterns):
    entries = []
    t1 = tk.Frame(root)

    for p in patterns:
        row = tk.Frame(root)
        names = ['Event1','Event2','Event3','Colors','Distance']
        for name in names:
            lab = tk.Label(row, text=name,width=10)
            lab.grid(row=0, column=names.index(name)+1)
        lab = tk.Label(row, width=20, text=p)
        lab.grid(row=1,column=0)
        number_of_events = 3
        for i in range(0,number_of_events):
            string_var = tk.StringVar(t1)
            if i == 0 and patterns.index(p) == 0:
                string_var.set('S')
            else:
                string_var.set('None')
            EV = ['S','D','HT','None']
            EV = tk.OptionMenu(row, string_var, *EV)
            EV.grid(row=1, column=i+1)
            entries.append((p+'_event'+str(i),string_var))
        colors = ['red','black','None']
        string_var = tk.StringVar(t1)
        string_var.set('None')
        entries.append((p + 'colors', string_var))
        colors_op = tk.OptionMenu(row, string_var, *colors)
        colors_op.config(width=10)
        colors_op.grid(row=1, column=number_of_events+1)
        string_var = tk.BooleanVar(t1)
        c = tk.Checkbutton(row,variable=string_var)
        entries.append((p + '_dis', string_var))
        c.grid(row=1, column=number_of_events+2)
        row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
    return entries

def floatValidation(S):
    flag = True
    for s in S:
        if s not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '.']:
            flag = False
    if flag:
        return True
    t3.bell()  # .bell() plays that ding sound telling you there was invalid input
    return False


def make_figure_parameters(self,root):
    bool_var = []
    c = []
    entries = []
    t1 = tk.Frame(root)
    names = ['color','labels','draw marked']
    for name in names:
        lab = tk.Label(t1, text=name,width=10)
        lab.grid(row=0, column=names.index(name))
        bool_var.append(tk.BooleanVar(t1))
        c.append(tk.Checkbutton(t1, variable=bool_var[names.index(name)],command=lambda param=bool_var[names.index(name)]: self.callback_color(param)))
        entries.append((name, bool_var[names.index(name)]))
        c[names.index(name)].grid(row=1, column=names.index(name))
        t1.pack(side=tk.TOP, padx=5, pady=5)

    row = tk.Frame(root)
    floatVal = (row.register(floatValidation), '%S')

    lab = tk.Label(row, width=2, text='x')
    ent = tk.Entry(row, validate='key', vcmd=floatVal,width=5)
    ent.insert("end",'10')
    lab.grid(row=0, column=0)
    ent.grid(row=0, column=1)
    entries.append(('x', ent))
    row.pack(side=tk.TOP, padx=0, pady=5)

    lab = tk.Label(row, width=5, text='X')
    lab.grid(row=0, column=2)

    lab = tk.Label(row, width=2, text='y')
    ent = tk.Entry(row, validate='key', vcmd=floatVal,width=5)
    ent.insert("end",'20')
    lab.grid(row=0, column=3)
    ent.grid(row=0, column=4)
    entries.append(('y', ent))
    row.pack(side=tk.TOP, padx=0, pady=5)
    return entries


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

    entries = []
    t3 = tk.Frame(root)
    intVal = (t3.register(intValidation), '%S')
    floatVal = (t3.register(floatValidation), '%S')
    lab = tk.Label(root, text='RSAM-finder')
    lab.config(width='50')
    lab.pack()
    for field in fields:
        if len(field) == 3:
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
            lab.grid(row=fields.index(field),column=0)
            ent.grid(row=fields.index(field),column=1)
            entries.append((field[0], ent))
            row.pack(side=tk.TOP, padx=0, pady=5)

    if 'general_settings' in fields:
        row = tk.Frame(root,background='grey')
        lab = tk.Label(row, text='create sigma from FASTA file', width=30,background='grey')
        lab.grid(row=0, column=0)
        create_sigma = tk.BooleanVar(t3)
        c = tk.Checkbutton(row, variable=create_sigma,background='black')
        c.grid(row=0, column=1)
        entries.append(('create_sigma', create_sigma))

        lab = tk.Label(row, text='Draw the gene tree', width=20,background='grey')
        lab.grid(row=0, column=2)
        draw = tk.BooleanVar(t3)
        c = tk.Checkbutton(row, variable=draw,background='black')
        c.grid(row=0, column=3)
        entries.append(('draw', draw))
        row.pack(side=tk.TOP, fill=tk.X, padx=0, pady=5)

        row = tk.Frame(root)
        lab = tk.Label(row, text='track solution number:', width=40)
        lab.grid(row=0, column=0)

        intVal = (t3.register(intValidation), '%S')
        ent = tk.Entry(row, validate='key', vcmd=intVal)
        entries.append(('track_solution', ent))
        ent.grid(row=0, column=1)
        row.pack(side=tk.TOP, fill=tk.X, padx=0, pady=5)
    return entries


class Main_Frame(object):
    def __init__(self, top=None):
        self.top = top
        img = ImageTk.PhotoImage(Image.open("logo.jpg"))
        panel = tk.Label(root, image=img)
        panel.pack(side="top", fill="both", expand="yes")

        parameter_lables = [['Specie tree extension', 'string', None], ['k', 'int', '50'],
                  ['Threshold Edges\nin Subtree', 'float', '0.1'],
                  ['HT event cost', 'float', '1'], ['Duplication event Cost', 'float', '1'],
                  ['Speciation event Cost', 'float', '0'], ['Loss event Cost', 'float', '0'], ['Gamma', 'float', '1'],
                  ['p', 'float', '0.05'], ['Number of Vertices to find', 'int', 1],'general_settings']

        nb = Notebook(self.top,'RSAM-Finder')
        t1 = nb.add_tab('General Prarmeters', parameters_frame, parameter_lables,'general')
        t2 = nb.add_tab('Pattern Definitions', parameters_frame,[],'patterns')
        t3 = nb.add_tab('Figure Prarmeters', parameters_frame, [], 'figure_params')

        self.b1 = tk.Button(top, text='Run',
                            command=lambda: self.bar_init(t1,t2,t3))
        self.b1.pack(side=tk.TOP, padx=5, pady=5)
        b2 = tk.Button(top, text='Quit', command=top.quit)
        b2.pack(side=tk.TOP, padx=5, pady=5)

        # run mainloop
        self.top.mainloop()

    def bar_init(self,t1,t2,t3):
        self.window = tk.Toplevel(root)
        ent1 = t1.work_task()
        ent2 = t2.work_task()
        ent3 = t3.work_task()
        msg = tk.Message(self.window, text='RSAM-finder in progress...')
        msg.grid(row=0, column=0)
        self.load_bar = ttk.Progressbar(self.window)
        self.load_bar.grid(row=1, column=0)
        # first layer of isolation, note var being passed along to the self.start_bar function
        # target is the function being started on a new thread, so the "bar handler" thread
        self.start_bar_thread = threading.Thread(target=self.start_bar, args=(ent1+ent2+ent3,))
        # start the bar handling thread
        self.start_bar_thread.start()

    def start_bar(self, ent):
        # the load_bar needs to be configured for indeterminate amount of bouncing
        self.load_bar.config(mode='determinate', maximum=100, value=0)
        # 8 here is for speed of bounce
        self.load_bar.start(8)
        # start the work-intensive thread, again a var can be passed in here too if desired
        self.work_thread = threading.Thread(target=self.work_task, args=(ent, ['p1','p2']))
        self.work_thread.start()
        # close the work thread
        self.work_thread.join()
        # stop the indeterminate bouncing
        self.load_bar.stop()
        # reconfigure the bar so it appears reset
        self.load_bar.config(value=0, maximum=0)
        self.window.destroy()
        msg = tk.Label(self.top, width=30, text='Results can be found in /data/results.')
        msg.pack()
        self.b1.destroy()

    def work_task(self, entries, patterns):
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
        if not p2_EV is None:
            p2 = (p2_EV, p2_colors, p2_dis)
        else:
            p2 = (None, None, False)
        vars = []
        for entry in entries:
            if entry[0].find('p1') == -1 and entry[0].find('p2') == -1:
                text = entry[1].get()
                vars.append(text)
        speciesTreespecification, k, TH_edges, HT_cost, D_cost, S_cost, loss_cost, gamma, p, number_of_planted_vertices, create_sigma, draw, track_solution,color,labels,draw_marked,x,y = vars
        if track_solution == "":
            track_solution = False
        RSAMfinder.main(speciesTreespecification, int(k), Decimal(TH_edges), int(HT_cost), int(D_cost), int(S_cost),
                        int(loss_cost),
                        Decimal(gamma), Decimal(p), int(number_of_planted_vertices), p1, p2, True, create_sigma,
                        track_solution, draw, color,labels,draw_marked,x,y)


class Notebook:
    def __init__(self, top,title):
        self.root = top
        self.root.title(title)
        self.notebook = ttk.Notebook(self.root)

    def add_tab(self, title, func,labels,specification):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text=title)
        tab = func(frame,labels,specification)
        self.notebook.pack()
        return tab

    def run(self):
        self.root.mainloop()


root = tk.Tk()
Main_Frame(root)


