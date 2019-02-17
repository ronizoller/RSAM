from datetime import datetime
import utiles

now = datetime.now()
for i in range(0,10000000):
    print('wait')
print(utiles.avg_dict_datetime({1:datetime.now()-now}))

