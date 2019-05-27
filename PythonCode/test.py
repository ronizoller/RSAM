import tkinter as tk
import RSAMfinder

patterns = ['Pattern 1: ', 'Pattern 2: ']

def run(entries):
    for e in entries:
        print('%s %s' %(str(e[0]),str(e[1].get())))
    quit()
    vars = []
    for entry in entries:
        text  = entry[1].get()
        vars.append(text)
    speciesTreespecification, k, TH_edges, HT_cost, D_cost, S_cost, gamma, p, number_of_planted_vertices, p1, p2 = vars
    RSAMfinder(speciesTreespecification,k,TH_edges,HT_cost,D_cost,S_cost,gamma, p,number_of_planted_vertices,  p1, p2,True)

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
        ind = patterns.index(p)
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
        c.grid(row=1,column=number_of_events+2)
        row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
    return entries

if __name__ == '__main__':
    root = tk.Tk()
    labels = [['Specie tree extension\n(will be used also for the S_edgelist extenstion)','string',None],['k','int','50'],['Threshold Edges\nin Subtree','float','0.1'],
              ['HT event cost','float','1'],['Duplication event Cost','float','1'],
              ['Speciation event Cost','float','0'],['Gamma','float','1'],
          ['p','float','0.05'],['Number of Vertices to find','int',None]]

    ents = makeform(root, labels)
    v = tk.BooleanVar()
    MODES = [("Double-Mode", True),("Single-Mod", False)]
    for text, mode in MODES:
        b = tk.Radiobutton(root, indicatoron=0, text=text, variable=v, value=mode)
        b.pack(side=tk.LEFT, padx=5, pady=5)

    b1 = tk.Button(root, text='Run',
                  command=(lambda e=ents: run(e)))
    b1.pack(side=tk.LEFT, padx=5, pady=5)
    b2 = tk.Button(root, text='Quit', command=root.quit)
    b2.pack(side=tk.LEFT, padx=5, pady=5)


    root.mainloop()
