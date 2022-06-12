

from path import Path

import glob



workdir = "/mnt/home/huanggy/project/20220325"

tmp = Path(workdir).glob('result/*/kraken2/*.report.mpa.txt')
print("tmp:===============>",tmp)
for f in tmp:
    # print(dir(f))
    print(f.basename(),type(f.basename()))

    if f.basename().find("21JS944002") != -1:
        print("find ======>yinxingkongzhi:",f.basename())

#print(tmp.__next__())
#print(tmp.__next__())


# a = glob.glob(Path(workdir)/'result/*/kraken2/*.report.mpa.txt')
# print(a)
# for f in a :
#     print(f)