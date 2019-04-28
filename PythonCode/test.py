def index_with_repeting(lst,item,list_of_returned):
    lst.reverse()
    index = lst.index(item)
    if list_of_returned[index] == 1:
        while index < len(list_of_returned) and list_of_returned[index] == 1:
            index += 1
    list_of_returned[index] = 1
    lst.reverse()
    return index,list_of_returned


number_of_planted_vertices = 10
marked_nodes = {'u597': [1,1, 'Double-mode'],'u5d97': [1,1, 'Double-mode']}
p1 = [1,1]
p2 = [1,1]

max_list = sorted([f + g for f in p1 for g in p2])
list_of_returns = [0]*len(max_list)

for nd,x in marked_nodes.items():
    itm = x[0] + x[1]
    ind, list_of_returns = index_with_repeting(max_list, itm, list_of_returns)
    print('index = '+ str(ind))

max_list.reverse()
for c, value in enumerate(max_list, 0):
    print(c, value,end=', ')

