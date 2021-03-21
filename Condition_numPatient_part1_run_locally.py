#!/usr/bin/env python
# coding: utf-8

# In[4]:


import pandas as pd
import numpy as np
import sys
import gc
import time
import itertools
import functools
import operator
import json
import psycopg2
import pickle
import networkx as nx
from networkx.readwrite import json_graph
import os

pd.set_option('display.max_colwidth', None)


# In[5]:


# geneate directories
parent_dir = "test/"
dirs = ["DG_for_colab", "DG_results", 
         "inclusions_with_patient", "inclusions_without_patient", 
         "exclusions_with_patient", "exclusions_without_patient",
        "frequency_tables", "frequency_tables_expand_by_f_and_num_criteria",
        "frequency_tables_expand_sorted_by_f", "condition_with_no_trials_and_patterns" ]
for dir in dirs:
    try:
        os.mkdir(parent_dir + dir)
    except OSError as error:  
        print(error)   


# In[6]:


conn = psycopg2.connect(
    host="localhost",
    database="cdm_5.2.2_SynPUF5pct",
    user="postgres",
    password="123456")
cur = conn.cursor()


# In[7]:


def compute_total_numPatients_in_database():
    sql = "select count(distinct person_id) from person"
    cur.execute(sql)
    result = cur.fetchone()
    if result != None:
        return result[0]
    else:
        return 0
compute_total_numPatients_in_database()


# In[8]:


## Setting
TOP_K_CRITERIA = 25
UPPER_BOUND = 5 #The upper bound of the number of criteria in a pattern
TOP_M_PATTERNS = 20
condition_ids = [313217]#, 317576, 314658, 314665] #Conditions that will be analyzed
condition_names = {313217:'Atrial fibrillation', 
                   317576:'Coronary arteriosclerosis',
                   314658:'Cardiomegaly',
                   314665:'Atrial flutter'}
general_criteria = [4274025, 4322976] #Some general criteria that will be excluded from the descendant_concept list.
TOTAL_NUMPATIENTS = compute_total_numPatients_in_database() #We may need to round it
color_dic = {1:'lightcoral', 2:'moccasin', 3:'palegreen', 4:'lightskyblue', 5:'violet'} #color of the nodes in the directed graph


# In[9]:


dic_domain_table = {"Condition": ("condition_occurrence", "condition_concept_id"), 
       "Device": ("device_exposure", "device_concept_id"), 
       "Drug": ("drug_exposure", "drug_concept_id"),
      "Observation": ("observation", "observation_concept_id"),
      "Measurement": ("measurement","measurement_concept_id"),
      "Procedure": ("procedure_occurrence", "procedure_concept_id")}


# In[10]:


def save_obj(obj, name):
    with open(''+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open('' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


# In[11]:


def concept_set_sql_formulate(concept_id, inc):
    #This function is to formulate the sql query to get all the peron_ids corresponding to a concept and its 
    #descendant concepts. The person in the list has at least one record of the concept or its descndant concepts.
    #It return the sql query and the concept_name of the concept.
    
    #Retrieve all the descendant_concept_ids of a concept
    sql1 = "SELECT descendant_concept_id     FROM concept_ancestor     where ancestor_concept_id = "+str(concept_id)+";"
    cur.execute(sql1)
    desc_ids = cur.fetchall()
    desc_ids = tuple([item[0] for item in desc_ids])
    
    
    print("number of descendants:", len(desc_ids))
    if len(desc_ids)>10000:
        print(concept_id, "number of descendant_concepts is larger than 10000")
    sql2 = ""
    dic_domain_concept = {"Condition":[], "Device":[],"Drug":[], "Observation":[], "Measurement":[],"Procedure":[]}
    concept_name = ""
    for desc_id in desc_ids:
        #Get the domain of the descendant_concepts, and categorize the descendant_cocnepts with their domain type.
        cur.execute(
                       "SELECT domain_id, concept_name FROM concept WHERE concept_id = "+str(desc_id)+";"
                        )
        result = cur.fetchone()
        if result != None:
            #If there is record of this descendant_concept in the database:
            domain = result[0]
            dic_domain_concept[domain].append(desc_id)
            if desc_id == concept_id:
                concept_name = result[1]
                print("concept_name:",concept_name)
    
    flag_domain_voc_contained = False
    for domain in dic_domain_concept.keys():
        table_name = dic_domain_table[domain][0]
        concept_column = dic_domain_table[domain][1]
        if dic_domain_concept[domain]!=[]:
            #If there exists descendant_concepts with the specific domain type, we query the table of this domain
            #and get all the person_ids with one of these descendant_concepts. 
            flag_domain_voc_contained = True
            sql2 += "(SELECT DISTINCT person_id FROM "+ table_name +                        " WHERE "+concept_column +" in ("
            id_list = functools.reduce(lambda x, y:str(x) +", "+str(y),dic_domain_concept[domain])
            sql2 += str(id_list)
            sql2 += ")) UNION "

   
    if flag_domain_voc_contained:
        sql2 = sql2[:-6]
           
           
    if inc == False:
        #If the concept is an exclusion:
        sql2 = "SELECT DISTINCT person_id FROM person EXCEPT ("+sql2+")"       
                
    return sql2, concept_name

def search_database_concept_set_descendant(concept_id, inc):
    #This function is to retrieve all the set of person_ids related to a concept or its descendant cocnepts. 
    #The person in the list has at least one record of the concept or its descndant concepts.
    sql, concept_name = concept_set_sql_formulate(concept_id, inc)
    if sql != "":
        cur.execute(sql)
        all_person_ids = cur.fetchall()
        all_person_ids = set(person_id[0] for person_id in all_person_ids)
    else:
        all_person_ids = set()
    return all_person_ids, concept_name
        


# In[12]:


def name_combination(row):
    ids = row['pattern']
    names = []
    for id in ids:
        if id>0:
            name = dic_inc_concept_names[id]
            names.append(name)
        else:
            name = dic_exc_concept_names[abs(id)]
            names.append("-"+name)
    return names
    

def flag_and_negate_exc(row):
    #This function is to flag whether the criterion is in the Top K inclusion or exclusion criteria list. 
    #It also flags the exclusion criterion.
    concept_id = row['criteria_concept_id']
    if int(concept_id in inclusions):
        row['flag']=1
    elif int(concept_id in exclusions):
        row['flag']=2
        #change the exclusion_criteria_id from positive to negative so as to distinguish exclusion and inclusion criteria. 
        row['criteria_concept_id'] = -row['criteria_concept_id']
    else:
        row['flag']=0
    return row
    


# In[13]:


def traverse_combinations(row):
    #This function is to traverse all the combinations of the criteria in a pattern, and assign them the frequency
    #with value equal to the one of the pattern.
    def position_id_conversion(comb):
        return tuple(map(lambda x: pattern[x], comb))
    
    pattern = eval(row["pattern"])
    num_criteria = len(pattern)
    id_combs = []
    upper_bound = min(UPPER_BOUND, num_criteria) #traverse the combinations of variables, the number of which is no more than UPPER_BOUND
    for k in np.arange(1, upper_bound+1):
        combinations = list(itertools.combinations(range(num_criteria),k))
        id_combs_k = map(position_id_conversion, combinations)
        id_combs += id_combs_k
    num_combs = len(id_combs)
    frequency = num_combs * [row[1]]
    return list(zip(id_combs, frequency))


# In[14]:


def compute_count(concept_ids):
    #Compute patient count of patterns and sub-patterns
    inter_persons = criteria[abs(concept_ids[0])]
    if len(concept_ids)>1:
        for id in concept_ids[1:]:
            inter_persons = inter_persons.intersection(criteria[abs(id)])
    count = len(inter_persons)
    return count
        
def add_sub_pattern_with_largest_reduction_rate(pattern): 
    #Draw the directed graph to represent the network topology of patterns and sub-patterns with the largest reduction rate.
    if len(pattern) > 1:
        #print(pattern)
        sub_patterns = list(itertools.combinations(pattern,len(pattern)-1))
        row = []
        count = 0 
        i = 0
        #print("length sub patterns", len(sub_patterns))
        pattern_with_largest_reduction_rate = ()
        while True:
            sub_pattern = sub_patterns[i]
            row_temp = freq_table.loc[freq_table['pattern']==str(sub_pattern)]
            count_temp = compute_count(sub_pattern)
            if i == 1:
                row = row_temp
                count = count_temp
                pattern_with_largest_reduction_rate = sub_pattern
            else:
                if count_temp>count:
                    row = row_temp
                    count = count_temp
                    pattern_with_largest_reduction_rate = sub_pattern

            i = i + 1
            if(i >= len(list(sub_patterns))):
                break
       # if str(pattern_with_largest_reduction_rate) in top_patterns:
       #     color = 'red'
       # else:
       #     color = 'blue'
        num_criteria = len(pattern_with_largest_reduction_rate)
        DG.add_nodes_from([(pattern_with_largest_reduction_rate,
                           {"names": row['pattern_concept_names'].values[0], 
                            "frequency": row['frequency'].values[0], 
                            "relative frequency": row['relative_frequency'].values[0],
                            "count": count,
                           "color": color_dic[num_criteria]})])
        count_pattern = DG.nodes[pattern]['count']
        reduction_rate = round((count - count_pattern)/count,5)
        reduction_rate = "{:.2%}".format(reduction_rate)
        #print(pattern, pattern_with_largest_reduction_rate, reduction_rate)
        DG.add_weighted_edges_from([(pattern_with_largest_reduction_rate, pattern, reduction_rate)])
        add_sub_pattern_with_largest_reduction_rate(pattern_with_largest_reduction_rate)

            
            
def add_sub_pattern(pattern): 
    #Draw the directed graph to represent the network topology of patterns and all sub-patterns.
    row = freq_table.loc[freq_table['pattern']==str(pattern)]
    count = compute_count(pattern)
    #if str(pattern) in top_patterns:
    #    color = 'red'
    #else:
    #    color = 'blue'
    num_criteria = len(pattern)
    DG.add_nodes_from([(pattern,
                           {"names": row['pattern_concept_names'].values[0], 
                            "frequency": row['frequency'].values[0], 
                            "relative frequency": row['relative_frequency'].values[0],
                            "count": count,
                           "color": color_dic[num_criteria]})])
    if len(pattern) != 1:
        for sub_pattern in itertools.combinations(pattern,len(pattern)-1):
            add_sub_pattern(sub_pattern)
            count_pattern = DG.nodes[pattern]['count']
            count_sub_pattern = DG.nodes[sub_pattern]['count']
            reduction_rate = round((count_sub_pattern - count_pattern)/count_sub_pattern,5)
            print(pattern, sub_pattern, reduction_rate)
            DG.add_weighted_edges_from([(sub_pattern, pattern, reduction_rate)])
        
    


# In[15]:


def get_table_with_numPatients_for_each_concept(domain):
    #Get the number of patients for each concept with domain, and create the table with field "concept_name", 
    #"concept_id" and "number of patients".
    table_name = dic_domain_table[domain][0]
    concept_id_column = dic_domain_table[domain][1]
    sql = "select p.concept_name, p.concept_id, count(p.person_id) as num_patients     from     (select distinct concept_name, concept_id, person_id     FROM "+ table_name +" JOIN concept     on concept_id = "+concept_id_column +" ) as p     group by p.concept_id, p.concept_name     order by count(p.person_id) desc;"
    cur.execute(sql)
    table_content = cur.fetchall()
    table_concept_numPatients = pd.DataFrame().from_dict(table_content)
    table_concept_numPatients.columns = ["concept_name", "concept_id", "num_patients"]
    return table_concept_numPatients


# In[16]:


"""
#Get the number of patients for each condition
condition_numPatients = get_table_with_numPatients_for_each_concept("Condition")
condition_numPatients.to_csv("test/condition_numPatients.csv", index=False, sep=",")

#condition_numPatient = pd.read_csv("test/condition_numPatients.csv",sep=",")
"""


# In[17]:


#Get the top K inclusion criteria and exclusion criteria. 
#If the number of criteria related to a condition is less than K, we keep all its criteria.
#If a criterion appears in both the top K inclusion and exclusion criteria, the exclusion criterion is removed. 
ctkb_all_disease_top_inclusion_criteria_50 = pd.read_csv("ctkb_all_disease_top_inclusion_criteria_50_1.csv",
                                                         sep=",",header=None)
ctkb_all_disease_top_exclusion_criteria_50 = pd.read_csv("ctkb_all_disease_top_exclusion_criteria_50_1.csv",
                                                         sep=",",header=None)
inc_top = {}
exc_top = {}
for i in range(len(ctkb_all_disease_top_inclusion_criteria_50)):
  inc = list(ctkb_all_disease_top_inclusion_criteria_50.iloc[i, 1:TOP_K_CRITERIA+1])
  inc = [int(item) for item in inc if ~np.isnan(item)]

  exc = list(ctkb_all_disease_top_exclusion_criteria_50.iloc[i, 1:TOP_K_CRITERIA+1])
  exc = [int(item) for item in exc if ~np.isnan(item)]

  inc_top[ctkb_all_disease_top_inclusion_criteria_50.iloc[i, 0]] = inc

  inter_set = set(inc).intersection(set(exc))
    
  exc_top[ctkb_all_disease_top_exclusion_criteria_50.iloc[i, 0]] = list(set(exc)-inter_set)


# In[18]:


#Preprocess the ctkb_all_criteria table
ctkb_all_criteria = pd.read_csv("ctkb_all_criteria.csv" )
ctkb_all_criteria["Count"] = 1
ctkb_all_criteria_dropped = ctkb_all_criteria.drop_duplicates()
ctkb_all_criteria_dropped = ctkb_all_criteria_dropped[ctkb_all_criteria_dropped.criteria_concept_id!="unmapped"]
ctkb_all_criteria_dropped.criteria_concept_id = pd.to_numeric(ctkb_all_criteria_dropped.criteria_concept_id)
#ctkb_all_criteria_dropped.head()

#Preprocess the ctkb_all_trials table
ctkb_all_trials = pd.read_csv("ctkb_all_trials.csv" )
ctkb_all_trials = ctkb_all_trials[["nctid","condition_concept_id"]]
ctkb_all_trials = ctkb_all_trials.drop_duplicates() 
#ctkb_all_trials.sort_values(by="nctid").head()


# In[19]:


#Compute the number of trials for each condition
dic_num_trials = ctkb_all_trials.groupby("condition_concept_id").count().to_dict()['nctid']
dic_str = json.dumps(dic_num_trials)
with open("test/num_trials_per_condition.txt", 'w') as writer:
    writer.write(dic_str)

"""
def jsonKeys2int(x):
    if isinstance(x, dict):
            return {int(k):v for k,v in x.items()}
    return x

with open("test/num_trials_per_condition.txt", 'r') as reader:
    dic=reader.read()
dic_num_trials = json.loads(dic,object_hook=jsonKeys2int)
"""


# In[20]:


for condition_id in condition_ids:
    print("condition:", condition_id)
    start = time.time()
    #####step 1: get the person_id set for each TOP K inclusion and exclusion criteria
    inclusions_original = inc_top[condition_id]
    exclusions_original = exc_top[condition_id]

    inclusions = []
    exclusions = []

    inclusions_no_patient = []
    exclusions_no_patient = []
    
    dic_inc_concept_names = {}
    dic_exc_concept_names = {}

    #Create the person_id set for each criterion in the Top K inclusion criteria list
    print("inclusion:")
    i = 1
    criteria  = {}
    for concept_id in inclusions_original:
        print("concept_id:",concept_id)
        if concept_id not in general_criteria:
            person_ids, concept_name = search_database_concept_set_descendant(concept_id=concept_id, inc=True)
            #Get the person_id set related to a condition and its descendant concepts
            print("patient_counts:", len(person_ids))
            if len(person_ids) != 0:
                criteria[concept_id] = person_ids
                inclusions.append(concept_id)
                dic_inc_concept_names[concept_id] = concept_name
            else:
                inclusions_no_patient.append(concept_id)
        print(i)
        i += 1
    
    #Create the person_id set for each criterion in the Top K exclusion criteria list
    print("exclusion:")
    i = 1
    for concept_id in exclusions_original:
        print("concept_id:",concept_id)
        if concept_id not in general_criteria:
            person_ids, concept_name = search_database_concept_set_descendant(concept_id=concept_id, inc=False)
            #Get the person_id set related to a condition and its descendant concepts
            print("patient_counts:", len(person_ids))
            if len(person_ids) != TOTAL_NUMPATIENTS:
                exclusions.append(concept_id)
                criteria[concept_id] = person_ids
                dic_exc_concept_names[concept_id] = concept_name

            else:
                exclusions_no_patient.append(concept_id)
        print(i)
        i += 1

    #save_obj(criteria, "test/patient_count_per_concept/condition_"+str(condition_id))
    save_obj(inclusions, "test/inclusions_with_patient/condition_"+str(condition_id))
    save_obj(exclusions, "test/exclusions_with_patient/condition_"+str(condition_id))
    save_obj(inclusions_no_patient, "test/inclusions_without_patient/condition_"+str(condition_id))
    save_obj(exclusions_no_patient, "test/exclusions_without_patient/condition_"+str(condition_id))

    
    
    ######step 2: Create table with pattern frequency for each condition
    
    #Join ctkb_all_trials table and ctkb_all_criteira_dropped table to 
    merged_df = ctkb_all_trials.merge(ctkb_all_criteria_dropped)
    dic_condition = dict(tuple(merged_df.groupby("condition_concept_id")))
    df = dic_condition[condition_id]
    
    #change the exclusion criteria id from positive to negative, and flag exclusion criterion in Top K list as 2, 
    #inclusion criterion in Top K list as 1, and criterion not in the list as 0.
    df = df.apply(flag_and_negate_exc, axis=1)
    df = df.loc[df['flag']!=0]
    df = df[['nctid','criteria_concept_id', 'Count']]
    df = df.sort_values(['nctid', 'criteria_concept_id'])
    
    #spread the table
    df = df.groupby(['nctid','criteria_concept_id'])['Count'].sum().unstack().fillna(0).astype(int)

    #generate patterns
    c_ids = np.array(df.columns)
    df["pattern"] = df.apply(lambda x: tuple(c_ids[np.where(x)]), axis=1)
    #df.to_csv("test/condition_tables/condition_"+str(condition_id)+".csv",
    #            header=True, sep="\t")

    #generate patern frequency for one condition
    freq_table = pd.DataFrame({'pattern':df['pattern'].tolist()})
    freq_table['frequency'] = 1
        #print(freq_table)
    new_table = freq_table.groupby(['pattern'])['frequency'].sum()
    new_table = new_table.reset_index()
    new_table['pattern_concept_names'] = new_table.apply(name_combination, axis=1)
    new_table.to_csv("test/frequency_tables/frequency_condition_"+str(condition_id)+".csv",
            header=True, sep="\t")

    del df, freq_table, new_table, merged_df
    gc.collect()
    df=pd.DataFrame()
    freq_table=pd.DataFrame()
    new_table=pd.DataFrame()
    merged_df = pd.DataFrame()
    dic_condition.clear()  



    #####step 3: find all the sub-patterns for each patterns, recaluate the pattern frequency, and generate sorted trial_frequency table
    num_trials = dic_num_trials[condition_id]
    freq_table = pd.read_csv("test/frequency_tables/frequency_condition_"+str(condition_id)+".csv",sep="\t", index_col=0)
    if freq_table.shape[0]!=0: #If the condition has related trials:
        #traverse all combinations of criteria in a pattern to get all sub-patterns
        temp = freq_table.apply(traverse_combinations, axis = 1).to_list()
        frequency_condition_expand = pd.DataFrame(functools.reduce(operator.iconcat, temp, []), columns=['pattern', 'frequency'])
        
        #recalculate the frequency of a pattern 
        new_table = frequency_condition_expand.groupby(['pattern'])['frequency'].sum()
        new_table = new_table.reset_index()
        new_table['relative_frequency'] = round(new_table['frequency']/num_trials,4)
        new_table['#_criteria'] = new_table.apply(lambda x: len(x['pattern']), axis=1)
        new_table['pattern_concept_names'] = new_table.apply(name_combination, axis=1)
        
        #Sort the table by frequency in descending order
        new_table = new_table.sort_values(['frequency'], ascending=[False], ignore_index=True)
        new_table[['pattern', 'pattern_concept_names', 'frequency', 'relative_frequency', '#_criteria']].to_csv("test/frequency_tables_expand_sorted_by_f/frequency_condition_"+str(condition_id)+"_expand_f.csv",
                    header=True, sep="\t")
        
        #Sort the table by number of criteria in the pattern in ascending order and frequency in descending order
        new_table['rank_in_group'] = new_table.groupby(by = ["#_criteria"]).cumcount()
        new_table.index=pd.MultiIndex.from_arrays([list(new_table['#_criteria']),list(new_table['rank_in_group'])])
        new_table = new_table.sort_index()                                
        new_table = new_table[['pattern', 'pattern_concept_names', 'frequency', 'relative_frequency']]
        new_table.to_csv("test/frequency_tables_expand_by_f_and_num_criteria/frequency_condition_"+str(condition_id)+"_expand_f_nc.csv",
                         header=True, sep="\t")
    else:
        #record the empty table
        freq_table.to_csv("test/condition_with_no_trials_and_patterns/frequency_condition_"+str(condition_id)+".csv",
                          header=True, sep="\t")
        
    del freq_table, new_table, temp
    gc.collect()
    freq_table = pd.DataFrame()
    temp = pd.DataFrame()
    new_table = pd.DataFrame()
    
    
    
    #####step 4: Compute the patient counts for each pattern, and draw the directed graph
    freq_table = pd.read_csv("test/frequency_tables_expand_by_f_and_num_criteria/frequency_condition_"+str(condition_id)+"_expand_f_nc.csv",sep="\t", index_col=[0,1])

    DG = nx.DiGraph()
    DG.add_nodes_from([('root',{"names":condition_names[condition_id],"color": "lightgray", "count": TOTAL_NUMPATIENTS})])
    for i in np.arange(UPPER_BOUND,0, -1):
        if i==1: #If we want to include all patterns with one criterion,
            freq_table_sub = freq_table.loc[[i]]
        else: #If we want to include Top M patterns with one criterion,
            freq_table_sub = freq_table.loc[(i,0):(i,TOP_M_PATTERNS-1)]
        for index, row in freq_table_sub.iterrows():
            pattern = eval(row['pattern'])

          ## add sub-pattern with largest reduction rate
            count = compute_count(pattern)
            num_criteria = len(pattern)
            DG.add_nodes_from([(pattern, 
                              {"names": row['pattern_concept_names'], 
                              "frequency": row['frequency'], 
                              "relative frequency": row['relative_frequency'],
                              "count": count,
                            "color": color_dic[num_criteria]})])
            add_sub_pattern_with_largest_reduction_rate(pattern)
            if i==1:
                count = DG.nodes['root']['count']
                count_pattern = DG.nodes[pattern]['count']
                reduction_rate = round((count - count_pattern)/count,5)
                reduction_rate = "{:.2%}".format(reduction_rate)
                DG.add_weighted_edges_from([('root', pattern, reduction_rate)])

      ##     
      ## add all sub-patterns 
      #add_sub_pattern(pattern)
    DG_json = json_graph.node_link_data(DG)
    patterns = [str(node['id']) for node in DG_json['nodes']]
    counts = [node['count'] for node in DG_json['nodes']]
      
    node_table = pd.DataFrame({'pattern': patterns, 'patient_counts': counts}) 
    freq_table = pd.read_csv("test/frequency_tables_expand_sorted_by_f/frequency_condition_"+str(condition_id)+"_expand_f.csv",sep="\t", index_col=[0])
    node_table = node_table.merge(freq_table, on="pattern", how='left')
    node_table = node_table.rename({"frequency":"trial_frequency", "relative_frequency": "trial_relative_frequency"}, axis='columns')
    node_table.to_csv("test/DG_results/top_"+str(TOP_M_PATTERNS)+"_nodes_condition_"+str(condition_id)+".csv", index=False,sep="\t")    
    
    edge_table = pd.DataFrame.from_dict(DG_json['links'])
    edge_table.columns = ['reduction_rate','source_pattern', 'target_pattern']
    edge_table = edge_table[['source_pattern', 'target_pattern', 'reduction_rate']]
    edge_table.to_csv("test/DG_results/top_"+str(TOP_M_PATTERNS)+"_edges_condition_"+str(condition_id)+".csv", index=False,sep="\t")
    
    save_obj(DG, "test/DG_for_colab/DG_top_"+str(TOP_M_PATTERNS)+"_condition_"+str(condition_id))
    
    del freq_table, node_table, edge_table
    gc.collect()
    freq_table = pd.DataFrame()
    node_table = pd.DataFrame()
    edge_table = pd.DataFrame()
    criteria.clear()
    DG.clear()
    
    
    end = time.time()
    print("")
    print("time:", end-start)
    print("")


# In[21]:


cur.close()


# In[ ]:





# In[ ]:


"""
def pattern_sql_formulate(pattern):
    This function is to generate sql query to compute the patient count for each pattern
    if pattern =="root": 
        sql = "SELECT COUNT(DISTINCT person_id) FROM person;"
        return sql
    i = 1
    sql = "SELECT COUNT(*) FROM ("
    for concept_id in pattern:
        if concept_id > 0 :
            sql +="("
            s,_ = concept_set_sql_formulate(concept_id, True)
            sql += s
            sql +=") INTERSECT "
        else:
            sql +="("
            s,_ = concept_set_sql_formulate(concept_id, True)
            sql += s
            sql +=") INTERSECT "
        i += 1
      
    sql = sql[:-10]
    sql +=") AS A;"
    return sql
"""

