#!/usr/bin/env python
# coding: utf-8

# In[1]:


import wget, os, bz2
import multiprocessing as mp
AP_counts = [60, 126]
counts = sum(AP_counts)
AP_sign = ["", "2"]


# In[2]:


def create_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)


# In[3]:


folder = "ssds/"
create_folder(folder)


# In[4]:


AP_name = "sqlApogee{}Object_{:04d}.csv.bz2"

for ii in range(2):
    for i in range(1, AP_counts[ii]+1):
        name = AP_name.format(AP_sign[ii], i)
        url = "https://dr17.sdss.org/sas/dr17/casload/apCSV/target/" + name
        download = True
        file = folder + name
        if os.path.exists(file):
            download = False
        if download:
            wget.download(url, out=file)
        print("\nFinished downloading {}".format(file))


# In[5]:


config = {}
config_file = open("../config.txt", 'r')
for l in config_file:
    l = l.strip()
    if len(l) > 0 and l[0] != "#":
        key, value = [x.strip() for x in l.split("=")]
        config[key] = value

import psycopg
table_name = "ssds"

# Connect to an existing database
with psycopg.connect("dbname={} host={} port={} user={} password={}".format(config["database"], config["hostname"], config["port"], config["username"], config["password"])) as conn:
    with conn.cursor() as cur:
        cur.execute("""DROP TABLE IF EXISTS {}""".format(table_name))
        cur.execute("""
        CREATE TABLE IF NOT EXISTS {} (
            {} BIGINT,
            j DOUBLE PRECISION,
            h DOUBLE PRECISION,
            k DOUBLE PRECISION,
            tmass_prox DOUBLE PRECISION)
        """.format(table_name, config["id_column"]))
        conn.commit()


# In[6]:


global_size = mp.Value("i", 0)
global_index = mp.Value("i", 0)
def generating_data(global_size, global_index):
    while True:
        local_index = -1
        ap_version, ap_index = -1, -1
        with global_index.get_lock():
            local_index = global_index.value
            global_index.value += 1
        if local_index >= counts:
            break
        for count in AP_counts:
            ap_version += 1
            if local_index < count:
                ap_index = local_index + 1
            else:
                local_index -= count
        cols = ["j", "h", "k", "tmass_prox"]
        file_name = folder + AP_name.format(AP_sign[ap_version], ap_index)
        f = bz2.open(file_name, mode='rt')
        size = -1
        col_indices = [-1 for i in range(len(cols))]
        data = []
        for l in f:
            tokens = l.split(',')
            if size == -1:
                for i, col in enumerate(cols):
                    for j, token in enumerate(tokens):
                        if col == token:
                            col_indices[i] = j
            else:
                vals = [float(tokens[j]) for j in col_indices]
                if min(vals) >= 0:
                    data.append(vals)
            size += 1
        f.close()
        start_index = -1
        with global_size.get_lock():
            start_index = global_size.value+1
            global_size.value += size
            print("File:", file_name, "Start_index:", start_index, "Size:", size)
        with psycopg.connect("dbname={} host={} port={} user={} password={}".format(config["database"], config["hostname"], config["port"], config["username"], config["password"])) as conn:
            with conn.cursor() as cur:
                with cur.copy("COPY ssds ({}, {}) FROM STDIN".format(config["id_column"], ','.join(cols))) as copy:
                    for i, row in enumerate(data):
                        copy.write_row([start_index + i] + row)
            conn.commit()
        print("File:", file_name, "finished inserting", size, "rows")


# In[ ]:


ps = [mp.Process(target = generating_data, args = (global_size, global_index)) for i in range(int(config["physical_core"]))]
for p in ps:
    p.start()
for p in ps:
    p.join()


# In[ ]:


print("Starts indexing ssds")
with psycopg.connect("dbname={} host={} port={} user={} password={}".format(config["database"], config["hostname"], config["port"], config["username"], config["password"])) as conn:
    with conn.cursor() as cur:
        cur.execute("""ALTER TABLE {} ADD PRIMARY KEY ({})""".format(table_name, config["id_column"]))
        conn.commit()
print("Finished indexing ssds")


# In[ ]:




