import math
import tree_operations_v1 as tree_operations
import utiles
import EfficiantVersion as effi
import os

def check_precentage(num,k):
    num = num/k
    if num >= 0 and num < 0.3:
        return 'red'
    elif num >= 0.3 and num < 0.6:
        return 'yellow'
    else : return 'green'

def first_n_digits(num, n):
    string_num = str(num)
    if string_num.find('e') == -1:
        return float(str(num)[:n+2])
    else:
        return 0

def fact(n, fact_lookup_table):
    f = math.factorial
    if n in fact_lookup_table:
        fact_n = fact_lookup_table[n]
    else:
        fact_n = f(n)
        fact_lookup_table.update({n: fact_n})
    return [fact_n,fact_lookup_table]

def normlize(l,largest):
    max = 0
    res = []

    for i in l:
        if i > max:
            max = i
    for i in l:
        if i == 1:
           res.append(1)
        else:
            res.append(i*largest/max)
    return res

def compute_avg_diff(l):
    sum = 0
    length = 0
    for i in l.values():
        if math.fabs(i[0] - i[1]) > 0:
            length += 1
        sum += math.fabs(i[0] - i[1])
    return sum/length

def p_value_calculation(fro, all, nCr_lookup_table, fact_lookup_table, accur, Pr_red, Pr_black, label):
    #print('Calculating p_value')
    if fro == all:
        return 0, nCr_lookup_table, fact_lookup_table
    else:
        res = 0
        for j in range(fro, all+1):
            temp, nCr_lookup_table, fact_lookup_table = utiles.nCr(all, j, nCr_lookup_table, fact_lookup_table)
            nCr_temp = utiles.first_n_digits(temp, accur)
            res = res + nCr_temp * pow(Pr_red, j) * pow(Pr_black, all - j)
        return res, nCr_lookup_table, fact_lookup_table

def init_internal_labels (tree, char, old_sigma, path):
    #print('Inisilasing internal leafs...')
    counter = 1
    dic = ''
    for nd in tree.postorder_node_iter():
        nd.label = char+str(counter)
        counter += 1
        if not nd.taxon == None:
            dic = dic + nd.label + ' : ' +  nd.taxon.label
            if char == 'u':
                dic = dic + ' (' + old_sigma[nd.taxon.label] + ') '
            dic = dic + '\n'
    #for nd in tree.postorder_node_iter():
    #   if nd.label == 'u365' or nd.label == 'u364':
    #       child = nd.child_nodes()
    #       nd.remove_child(child[0])
    #       print('     '+str(child[0])+' was deleted')
    to_create = path + '/saved_data/'
    os.makedirs(os.path.dirname(to_create), exist_ok=True)
    tree_operations.collapse_edges(tree)
    if char == 'u':
        file = open(path + '/saved_data/G_keys.txt', 'w')
        file.write(str(dic))
        file.close()
    if char == 'x':
        file = open(path + '/saved_data/S_keys.txt', 'w')
        file.write(str(dic))
        file.close()
    #print('Finished inisilasing internal leafs.\n')
    return tree

class heap_items:
    def __init__(self, cost,node1,node2):
        self.val = (cost,node1,node2)
    def __lt__(self, other):
        if (self.val[0]==other.val[0]):
            return True
        else:
            return self.val[0]<other.val[0]

def nCr(n,r, nCr_lookup_table, fact_lookup_table):
    #print('     Calculating '+str(n)+' choose '+str(r)+'...\n')
    if (n,r) in nCr_lookup_table:
        #print('     Finished calculating '+str(n)+' choose '+str(r))
        return [nCr_lookup_table[(n,r)],nCr_lookup_table,fact_lookup_table]
    else :
        temp = utiles.fact(n,fact_lookup_table)
        fact_n = temp[0]
        fact_lookup_table = temp[1]

        temp = utiles.fact(r, fact_lookup_table)
        fact_r = temp[0]
        fact_lookup_table = temp[1]

        temp = utiles.fact((n-r), fact_lookup_table)
        fact_nr = temp[0]
        fact_lookup_table = temp[1]

        nCr_lookup_table.update({(n,r): fact_n // fact_r // fact_nr})
        #print('Finished calculating '+str(n)+' choose '+str(r))
        return  nCr_lookup_table[(n,r)],nCr_lookup_table,fact_lookup_table

def expected(probs, l):
    zipped = list(zip(probs, l))
    res = 0
    for i in range(0, len(zipped)):
        res = res + zipped[i][0] * zipped[i][1]
    return res

def string_into_array (string):
    print ("[ '",end = "")
    for i in range (0,len(string)):
        if string[i] != " ":
            print(string[i],end = '')
        else:
            print("', '",end = "")
    print (" ]",end = "")
    quit()


def avg_dict_datetime (dict):
    sum_hours = 0
    sum_mini = 0
    sum_sec = 0
    length = 0
    for d,td in dict.items():
        days = td.days
        hours, remainder = divmod(td.seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        sum_hours += hours
        sum_mini += minutes
        sum_sec += seconds
        length += 1
    return round(sum_hours/length+(sum_mini/length)/60+(sum_sec/length)/6000,4)

def map_max(map,number_of_fields):
    max = 0
    for u,array in map.items():
            for k in range(0,2):
                if array[0] > max:
                    max = array[j][k]
    return max

def compare_dict_entries(dict):
    keys = []
    for k,internal_dict in dict.items():
        for u,couple in internal_dict.items():
            if u not in keys:
                keys.append(u)
    new_dic = {}
    for k,internal_dict in dict.items():
        new_internal_dic = {}
        for key in keys:
            if key in internal_dict:
                new_internal_dic.update({key:internal_dict[key]})
            else: new_internal_dic.update({key:(0,0)})
        new_dic.update({k:new_internal_dic})
    return new_dic

def flip_list(to_flip,length):
    res = []
    for i in range(0, length):
        res_temp = []
        for l in to_flip:
            res_temp.append(l[i])
        res = res + [res_temp]
    return res

def number_of_different_vertex_in_dict(dict):
    different = []
    res = 0
    for rand_num,list_of_scores in dict.items():
        for u,score in list_of_scores.items():
            if u not in different:
                different.append(u)
                res += 1
    return res


def average_of_list(dict,num):
    res = {}
    if dict == {}:
        return {}
    for rand_num, list_of_scores in dict.items():
        for u, score in list_of_scores.items():
            if not u in res:
                sum = [0,0]
                for rand_num1, list_of_scores1 in dict.items():
                    for u1, score1 in list_of_scores1.items():
                        if u1 == u:
                            sum[0] += score1[0]
                            sum[1] += score1[1]
                res.update({u: (sum[0]/num,sum[1]/num)})
    #print('         Average list: %s\n' % str(res))
    return res

def average_of_list_check_diff(dict):
    res = {}
    if dict == {}:
        return {}
    for TH, list_for_TH in dict.items():
        temp_res = []
        for list_for_iter in list_for_TH:
            for tup in list_for_iter:
                u = tup[0]
                if not u in temp_res:
                    temp_res.append(u)
        res.update({TH: temp_res})
    return res

def find_unmarked(all_marked_for_TH,G,RSAM):
    list_of_unmarked_TH = []
    if not RSAM:
        for v in G.postorder_node_iter():
            flag = False
            for u in all_marked_for_TH:
                if v.label == u[0]:
                    flag = True
            if not flag:
                list_of_unmarked_TH.append(v.label)
    else:
        for v in G.postorder_node_iter():
            flag = False
            for u,score in all_marked_for_TH.items():
                if v.label == u:
                    flag = True
            if not flag:
                list_of_unmarked_TH.append(v.label)
    return list_of_unmarked_TH

def calculate_presentage(all_marked_list,all_unmarked_list,planted_vertex):
    sensitivity = {}
    specifity = {}
    TN = {}
    for TH,unmarked_list in all_unmarked_list.items():
        TN.update({TH:[]})
        temp_TN = 0
        for u in unmarked_list:
            if u not in planted_vertex:
                temp_TN += 1
        TN[TH].append([temp_TN])
    for TH, list_for_TH_sol in all_marked_list.items():
        sensitivity.update({TH: []})
        specifity.update({TH: [] })
        i = 0
        planted_vertex_to_check = planted_vertex.copy()
        temp_FP = 0
        temp_TP = 0

        for u in list_for_TH_sol:
            u = u[0]
            if u not in planted_vertex:
                temp_FP += 1
            else:
                temp_TP += 1
                if u in planted_vertex_to_check:
                    planted_vertex_to_check.remove(u)
        temp_FN = len(planted_vertex_to_check)
        sensitivity[TH].append(temp_TP / (temp_TP + temp_FN))

        specifity[TH].append(TN[TH][i][0] / (TN[TH][i][0] + temp_FP))
        i += 1
    res = {}
    print('sen: '+str(sensitivity))
    for TH, sen_list in sensitivity.items():
        sensitivity_sum = 0
        for score in sen_list:
            sensitivity_sum += score
        res.update({TH:round(sensitivity_sum/len(sen_list),2)})
    for TH, spe_list in specifity.items():
        specifity_sum = 0
        for score in spe_list:
            specifity_sum += score
        res.update({TH: (1-specifity_sum / len(spe_list),res[TH])})
    return res

def frange(start,end,step):
    res = []
    i = start
    while i<= end:
        res.append(round(i,1))
        i += step
    return res

def kmin_positive (l,k,H,nodes_table):
    l_temp = l.copy()
    res = []
    for i in range(0,min([len(l_temp),k])):
        min_ind = [x[1]['cost'] for x in l_temp].index(min([x[1]['cost'] for x in l_temp]))
        temp_item = l_temp[min_ind]
        if temp_item[1]['cost'] == math.inf:
            return res
        l_temp.remove(l_temp[min_ind])
        if type(temp_item) == tuple:
            res.append(temp_item)
        else:
            res.append(
                effi.find_nodes_in_hypergraph(H, temp_item['s'], temp_item['t'], temp_item['list_place'], nodes_table)[0])
    #if len(l) < k:
        #res = res +[None]*(k-len(l))
    return res

def kmin_list(l_input,subtree1,subtree2,k,H,nodes_table):
    if l_input != []:
        node_index = l_input[0][0]
        l_new = (l_input[0][1]['l']).copy()
        l = [(node_index,item) for item in l_new] + subtree1.copy() + subtree2.copy()
        list_of_values = [nd[1]['cost'] for nd in l]
    else:
        l = subtree1.copy() + subtree2.copy()
        list_of_values = [x[1]['cost'] for x in subtree1 + subtree2]
    res = []
    #print('l: ' + str(l)+' k: '+str(k))
    end = len(list_of_values)
    for i in range(0,end):
        if min(list_of_values) == math.inf:
            return res
        else:
            temp_index = list_of_values.index(min(list_of_values))
            temp_item = l[temp_index]
            #print('temp_index: %s, min value: %s,list of values: %s\n\n' % (str(temp_index),str(temp_item),str(list_of_values)))
            l.remove(l[temp_index])
            list_of_values.remove(list_of_values[temp_index])
            if type(temp_item) == tuple:
                res.append(temp_item)
            else:
                res.append(
                    effi.find_nodes_in_hypergraph(H, temp_item['s'], temp_item['t'], temp_item['list_place'], nodes_table)[0])
    #if end < k:
        #res = res +[None]*(k-end)
    #print('min form l: '+str(res)+'\n')
    return res

def update_top_ranking_list(new_score,list):
    j = 0
    if new_score > list[0]:
        while new_score > list[j] and j < len(list)-1:
            j += 1
        if new_score < list[j]:
            j -= 1
        temp = list[:j+1]
        list[j] = new_score
        while j > 0:
            list[j-1] = temp[j]
            j -= 1
    return list

def index_with_repeting(lst,item,list_of_returned):
    lst.reverse()
    index = lst.index(item)
    if list_of_returned[index] == 1:
        while index < len(list_of_returned) and list_of_returned[index] == 1:
            index += 1
    if index == len(list_of_returned):
        return index -1,list_of_returned
    list_of_returned[index] = 1
    lst.reverse()
    return index,list_of_returned