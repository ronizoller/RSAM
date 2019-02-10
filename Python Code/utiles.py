import math
import tree_operations
import utiles

def check_precentage(num):
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
    print('Inisilasing internal leafs...')
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
    tree_operations.collapse_edges(tree)
    if char == 'u':
        file = open(path + '/saved_data/G_keys' + '.txt', 'w')
        file.write(str(dic))
        file.close()
    if char == 'x':
        file = open(path + '/saved_data/S_keys' + '.txt', 'w')
        file.write(str(dic))
        file.close()
    print('Finished inisilasing internal leafs.\n')
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

def string_into_array (string):
    print ("[ '",end = "")
    for i in range (0,len(string)):
        if string[i] != " ":
            print(string[i],end = '')
        else:
            print("', '",end = "")
    print (" ]",end = "")
    quit()

def map_max(map,number_of_fields):
    max = 0
    for u,array in map.items():
        for j in range(0,number_of_fields):
            for k in range(0,2):
                if array[j][k] > max:
                    max = array[j][k]
    print('map max = %s' % str(max))
    return max

def compare_dict_entries(dict):
    print(dict)
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
    print('         Average list: %s\n' % str(res))
    return res

def calculate_presentage(all_marked_list,planted_vertex):
    planted_vertex = planted_vertex[0]
    false_pos = {}
    false_neg = {}
    for index, marked_list in all_marked_list.items():
        planted_vertex_to_check = planted_vertex.copy()
        false_pos.update({index: 0})
        false_neg.update({index: 0})
        for u, score in marked_list.items():
            if u not in planted_vertex:
                false_pos[index] += 1
            if u in planted_vertex_to_check:
                planted_vertex_to_check.remove(u)
        if planted_vertex_to_check != []:
            false_neg[index] = len(planted_vertex_to_check)
        else:
            false_neg[index] = 0
        if len(marked_list) == 0:
            false_pos[index] = 0
        else:
            false_pos[index] = false_pos[index] / len(marked_list)
        false_neg[index] = false_neg[index] / len(planted_vertex)
    sum_false_pos = 0
    for i, score in false_pos.items():
        sum_false_pos += score
    sum_false_neg = 0
    for i, score in false_neg.items():
        sum_false_neg += score
    return [round((sum_false_pos / len(all_marked_list))*100,2), round((sum_false_neg / len(all_marked_list))*100,2)]
