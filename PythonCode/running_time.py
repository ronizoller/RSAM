import utiles
import Tree_Grnarator
from datetime import datetime
import os

times = {}
path = '/users/studs/bsc/2016/ronizo/PycharmProjects/RSAM/simulator_data/running_time'
for number_of_leafs in utiles.frange(100,1000,100):
    path_curr = path + '/number_of_leavs:'+number_of_leafs+'/'
    os.makedirs(os.path.dirname(path_curr), exist_ok=True)
    Tree_Grnarator.main(number_of_leafs,path_curr)
    times.update({number_of_leafs:[datetime.now(),0]})
