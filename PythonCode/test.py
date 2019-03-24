def update_top_ranking_list(new_score,list):
    print(new_score)
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
            print(list)
            print(temp)
            print(list[j-1])
            print('\n')
            j -= 1
    return list


top = [0] * 4

for num in [10,5,30,100,24444444,33,2134]:
    top = update_top_ranking_list(num,top)
print(top)