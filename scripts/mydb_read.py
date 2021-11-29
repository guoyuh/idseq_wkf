#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import pymysql

db = pymysql.connect(
    host="127.0.0.1",
    port=3306,
    user = "huanggy",
    password="654321",
    database = "taxdb",
    charset="utf8"
)

def read(tb_name):
    sql = "SELECT * FROM {0};".format(tb_name)
    rows, length = db.select(sql)
    data = []
    for row in rows:
        data.append(row)
    return data



read("userinfo")