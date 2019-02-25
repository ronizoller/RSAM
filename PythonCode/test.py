import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
import matplotlib.pyplot as plt

dict=  {0: [[('u259', [(0.0, 0.0007024793388429752), (3, 0), 'v blacks HT and reds under w'])], [('u259', [(0.0, 0.0007231404958677686), (3, 0), 'v blacks HT and reds under w'])], [('u259', [(0.0, 0.0007575757575757575), (3, 0), 'v blacks HT and reds under w'])], [('u259', [(0.0, 0.0007359307359307359), (3, 0), 'v blacks HT and reds under w'])], [('u259', [(0.0, 0.0007024793388429752), (3, 0), 'v blacks HT and reds under w'])]], 1: [[('u259', [(0.0, 0.0007024793388429752), (3, 0), 'v blacks HT and reds under w'])], [('u259', [(0.0, 0.0007231404958677686), (3, 0), 'v blacks HT and reds under w'])], [('u259', [(0.0, 0.0007575757575757575), (3, 0), 'v blacks HT and reds under w'])], [('u259', [(0.0, 0.0007359307359307359), (3, 0), 'v blacks HT and reds under w'])], [('u259', [(0.0, 7024793388429752), (3, 0), 'v blacks HT and reds under w'])]]}
res = {}
num = 5
for TH, list_for_TH in dict.items():
    temp_res = {}
    for list_for_iter in list_for_TH:
        for tup in list_for_iter:
            u = tup[0]
            if not u in temp_res:
                sum = [0, 0]
                for list_for_iter1 in list_for_TH:
                    for tup1 in list_for_iter1:
                        u1 = tup1[0]
                        if u1 == u:
                            sum[0] += tup1[1][0][0]
                            sum[1] += tup1[1][0][1]
                print(sum)
                temp_res.update({u: (sum[0] / num, sum[1] / num)})
    res.update({TH:temp_res})
print(res)