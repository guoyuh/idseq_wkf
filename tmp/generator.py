
import time

def test():
    d = {
    "a":1,
    "b":2,
    "c":3,
    "d":4,
    "e":5,
    "f":6,
    "g":7
    }
    for k,v in d.items():
        yield k



data = test()

for i in data:
    time.sleep(2)
    print(i)
    

    
