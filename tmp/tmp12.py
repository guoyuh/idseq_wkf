import json
import re
with open("/mnt/home/huanggy/project/dme00001.json","r") as f:
    fj = f.read()
    kojson = json.loads(fj)

with open("/mnt/home/huanggy/project/newKegg.tsv", "w") as k:
    for i in kojson['children']:
        
        ii = i['name'].replace(" ", "\t", 1)
        for j in i['children']:
            
            jj = j['name'].replace(" ", "\t", 1)
            for m in j['children']:
                
                if re.findall(r"ko\d{5}",m['name']):
                    mm = "ko" + m['name'].replace(" ", "\t", 1)
                else :
                    mm = m['name'].replace(" ", "\t", 1)
                try :
                    for n in m['children']:
                        
                        nn = n['name'].replace(" ", "\t", 1).replace("; ", "\t")
                        k.write(ii + "\t" + jj + "\t" + mm + "\t" + nn + "\n")
                except :
                    
                    nn = " \t \t "
                    k.write(ii + "\t" + jj + "\t" + mm + "\t" + nn + "\n")