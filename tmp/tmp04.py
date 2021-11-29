
import pandas as pd 
left = pd.DataFrame({'key1': ['K0', 'K0', 'K1', 'K2'],
                      'key2': ['K0', 'K1', 'K0', 'K1'],
                         'A': ['A0', 'A1', 'A2', 'A3'],
                         'B': ['B0', 'B1', 'B2', 'B3']})

right = pd.DataFrame({'key1': ['K0', 'K1', 'K1', 'K2'],
                      'key2': ['K0', 'K0', 'K0', 'K0'],
                         'C': ['C0', 'C1', 'C2', 'C3'],
                         'D': ['D0', 'D1', 'D2', 'D3']})

result = pd.merge(left, right, on=['key1', 'key2'])
# 同时传入两个Key，此时会进行以['key1','key2']列表的形式进行对应，left的keys列表是：[['K0', 'K0'],['K0', 'K1'],['K1', 'K0'],['K2', 'K1']],
#left的keys列表是：[['K0', 'K0'],['K1', 'K0'],['K1', 'K0'],['K2', 'K0']]，因此会有1个['K0', 'K0']、2个['K1', 'K0']对应。

print(left)
print(right)
#print(result)

result2 = pd.merge(left, right, how='outer', on=['key1', 'key2'])

print(result2)