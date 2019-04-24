def index_with_repeting(lst,item,list_of_returned):
    index = lst.index(item)
    if list_of_returned[index] == 1:
        while index < len(list_of_returned) and list_of_returned[index] == 1:
            index += 1
    list_of_returned[index] = 1
    return index,list_of_returned

list1 = ['a','a','b']
list_of_returned = [0]*len(list1)
for lett in ['b','a','a']:
    ind,list_of_returned = index_with_repeting(list1, lett, list_of_returned)
    print(list_of_returned)
    print(ind)