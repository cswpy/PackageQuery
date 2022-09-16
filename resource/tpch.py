#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os, subprocess
import multiprocessing as mp


# In[2]:


config = {}
config_file = open("../config.txt", 'r')
for l in config_file:
    l = l.strip()
    if len(l) > 0 and l[0] != "#":
        key, value = [x.strip() for x in l.split("=")]
        config[key] = value

import psycopg
table_name = "tpch"

# Connect to an existing database
with psycopg.connect("dbname={} host={} port={} user={} password={}".format(config["database"], config["hostname"], config["port"], config["username"], config["password"])) as conn:
    with conn.cursor() as cur:
        cur.execute("""DROP TABLE IF EXISTS {}""".format(table_name))
        cur.execute("""
        CREATE TABLE IF NOT EXISTS {} (
            {} BIGINT,
            quantity DOUBLE PRECISION,
            price DOUBLE PRECISION,
            discount DOUBLE PRECISION,
            tax DOUBLE PRECISION)
        """.format(table_name, config["id_column"]))
        conn.commit()


# In[3]:


if not os.path.exists("tpch-kit"):
    result = subprocess.run("""
        git clone https://github.com/gregrahn/tpch-kit.git;
        cd tpch-kit/dbgen;
        make MACHINE=LINUX DATABASE=POSTGRESQL
    """, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode == 0:
        print("git clone and make succesffully")
    else:
        print("Error when executing git clone and make")


# In[4]:


sf = 300
thread = 44
global_size = mp.Value("i", 0)
def generating_data(global_size, sf, thread):
    file_name = "tpch-kit/dbgen/lineitem.tbl.{}".format(thread)
    if not os.path.exists(file_name):
        result = subprocess.run("""
            cd tpch-kit/dbgen; 
            ./dbgen -s {} -S {} -C {} -T L
        """.format(sf, thread, sf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            print("Error when executing dbgen for thread", thread)
    if os.path.exists(file_name):
        size = 0
        data = []
        for row in open(file_name):
            q, p, d, t = [float(x) for x in row.split("|")[4:8]]
            d *= p
            t *= p
            data.append([q, p, d, t])
            size += 1
        start_index = -1
        with global_size.get_lock():
            start_index = global_size.value+1
            global_size.value += size
            print("Thread:", thread, "Start_index:", start_index, "Size:", size)
        with psycopg.connect("dbname={} host={} port={} user={} password={}".format(config["database"], config["hostname"], config["port"], config["username"], config["password"])) as conn:
            with conn.cursor() as cur:
                with cur.copy("COPY tpch ({}, quantity, price, discount, tax) FROM STDIN".format(config["id_column"])) as copy:
                    for i, row in enumerate(data):
                        copy.write_row([start_index + i] + row)
            conn.commit()
        print("Thread", thread, "finished inserting", size, "rows")


# In[10]:


size = int(config["main_memory_size"]) // 4
for i in range(1, sf, size):
    a = i
    b = min(i + size, sf)
    print("Start collecting range ({},{})".format(a, b))
    ps = [mp.Process(target = generating_data, args = (global_size, sf, i)) for i in range(a, b+1)]
    for p in ps:
        p.start()
    for p in ps:
        p.join()


# In[37]:


print("Starts indexing tpch")
with psycopg.connect("dbname={} host={} port={} user={} password={}".format(config["database"], config["hostname"], config["port"], config["username"], config["password"])) as conn:
    with conn.cursor() as cur:
        cur.execute("""ALTER TABLE {} ADD PRIMARY KEY ({})""".format(table_name, config["id_column"]))
        conn.commit()
print("Finished indexing tpch")


# In[ ]:




