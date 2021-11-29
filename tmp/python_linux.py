import subprocess

cmd = "cat /mnt/home/huanggy/project/20211029/result/filter/LX2110965.f.fastq | wc -l " 
ret = subprocess.run(["wc","-l","/mnt/home/huanggy/project/20211029/result/filter/LX2110965.f.fastq"])
print(ret)


print("2222222222222222222222222")
ret2 = subprocess.run(["wc","-l","/mnt/home/huanggy/project/20211029/result/filter/LX2110965.f.fastq"])
print(ret2)

print('333333333333333333333333333')
ret3 = subprocess.getstatusoutput(cmd)
print(ret3)

print('44444444444444444444444')
ret4 = subprocess.getoutput(cmd)
print(ret4)

print('5555555555555555555555555')
ret5 = subprocess.check_call(cmd,shell=True)
print(ret5)

print('666666666666666666666666')
ret6 = subprocess.check_output(cmd,shell = True)
print(ret6)



a = ' LX2112143FFPED102G5kxT1 '
print(a)
print(a.lstrip().rstrip())