import datetime
import numpy as np
import hafnian
import sys

def sqbs_arr(m):
    """
    return str of Beamsplitter arrangement in "1,4,7,12,14" format
    """
    barr = ""
    for j in range(1,m):
        bs_arrj = ""
        if j%2==1:
            for k in range(1,m,2):
                if j!=1 or k!=1:
                    bs_arrj = bs_arrj + ","
                bs_arrj = bs_arrj + str(k)
        else:
            for k in range(2,m-1,2):
                if j!=1:
                    bs_arrj = bs_arrj + ","
                bs_arrj = bs_arrj + str(k)
        barr = barr + bs_arrj
    return barr

def timestamp():
    """
    return timestamp string: yyyy-mm-dd-hh-mmss
        e.g: "2020-05-20-815" or "2020-05-20-19"
    """
    date = str(datetime.date.today())
    hr = str(datetime.datetime.now().hour)
    minu = str(datetime.datetime.now().minute)
    sec = str(datetime.datetime.now().second)
    time = date+"-"+hr+"-"+minu+""+sec #timestamp file
    return time

def save():
    pass

def load():
    pass


