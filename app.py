import streamlit as st
from threading import activeCount
import matplotlib.pyplot as plt
import pickle
import umap
import io
import numpy as np
import pandas as pd

dfcurrent=pd.read_csv("G25_Current_DNA.csv")
Xcurrent=dfcurrent.drop(columns=['DNA sample ethnicity and id','DNA sample ethnicity','sample id'])

dfancient=pd.read_csv("G25_Ancient_DNA.csv")
Xancient=dfancient.drop(columns=['DNA sample ethnicity and id','DNA sample ethnicity','sample id'])
#Xcurrent

c=pd.read_csv('clustergmm15.csv')
c=c.drop(columns=['Unnamed: 0'])

frames = [dfcurrent, dfancient]
dfcombined = pd.concat(frames)
Xcombined=dfcombined.drop(columns=['DNA sample ethnicity and id','DNA sample ethnicity','sample id'])
#dfcombined

dfcurrentgroup=dfcurrent.groupby(['DNA sample ethnicity']).mean().reset_index()
Xcurrentgroup=dfcurrentgroup.drop(columns=['DNA sample ethnicity'])

dfadnalineages=pd.read_csv("adnalineages.csv")
dfancientpcadna = pd.merge(dfadnalineages,dfancient)
dfancienthpg=dfancientpcadna.groupby(['Assigned Mutation']).mean().reset_index()


uploaded_file = st.file_uploader("Choose a file")
if uploaded_file is not None:
     # To read file as bytes:
     bytes_data = uploaded_file.getvalue()
     st.write(bytes_data)

     # To convert to a string based IO:
     stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
     st.write(stringio)

     # To read file as string:
     string_data = stringio.read()
     st.write(string_data)

     # Can be used wherever a "file-like" object is accepted:
     input = pd.read_csv(uploaded_file)
     st.write(input)
    
    
def euclidean_distance(p, q):
    n = 25 # dimensions
    return sum([(p[i] - q[i]) ** 2 for i in range(n)]) ** 0.5

p=Xcombined.iloc[735:740]

dfdistances=dfcombined
distances=[]
#Induvidual

for i in range(len(Xcombined)):
  distances.append(euclidean_distance(Xcombined.iloc[i],p.iloc[3]))
dfdistances['distances']=distances

dfdistances=dfdistances.sort_values(by=['distances'])

p=Xcombined.iloc[735:740]

import numpy as np
p1 = np.zeros((len(p),len(c)))
p2 = np.zeros((len(p),len(c)))

for i in range(len(p)):
  for j in range(len(c)):
    p1[i][j]=euclidean_distance(c.iloc[j],p.iloc[i])
distmat=pd.DataFrame(p1)

for q in range(len(p)):
  tot=distmat.iloc[q].sum()

  for w in range(len(c)):
    p2[q,w]=1-((tot-distmat.iloc[q,w])/tot)

ancestry=pd.DataFrame(p2)


filename = 'umap_model.sav'

#input=Xancient.iloc[2288:2316]
loaded_model = pickle.load(open(filename, 'rb'))
test_embedding = loaded_model.transform(input)
print(test_embedding)

dfucurrentgroup=pd.DataFrame(ucurrentgroup)
test_embedding = fit.transform(input)
dftest_embedding=pd.DataFrame(test_embedding)


#dfdistances
#print(dfdistances['DNA sample ethnicity and id'].iloc[:3])

Tools = st.selectbox("Choose your Tool", ["Distance Tool", "PCA(Principal Component Analysis) Tool","ML Ancestry Tool","Ancient DNA Lineage Tool","Umap Plot Tool"]) 

if Tools == "Distance Tool":
  st.dataframe(dfdistances)


elif Tools == "ML Ancestry Tool":
  st.dataframe(ancestry)


elif Tools == "PCA(Principal Component Analysis) Tool":
  fig, ax = plt.subplots(figsize=(45, 99))
  ax.scatter(dfcurrentgroup['1'], dfcurrentgroup['2'],s = 1)
  #ax.scatter(point['1'],point['2'],s=500)

  for i in range(len(dfcurrentgroup)):
    ax.annotate(dfcurrentgroup['DNA sample ethnicity'][i], (dfcurrentgroup['1'][i], dfcurrentgroup['2'][i]))
  st.pyplot()


elif Tools == "Ancient DNA Lineage Tool":

  fig2, ax2 = plt.subplots(figsize=(30, 24))
  ax2.scatter(dfancienthpg['1'], dfancienthpg['2'],s = 1)
  #ax.scatter(point['1'],point['2'],s=500)
  for i in range(len(dfancienthpg)):
    ax2.annotate(dfancienthpg['Assigned Mutation'][i], (dfancienthpg['1'][i], dfancienthpg['2'][i]))
  st.pyplot()


elif Tools == "Umap Plot Tool":
    fig3, ax3 = plt.subplots(figsize=(50, 200))
    scatter3 = ax3.scatter(dfucurrentgroup[0], dfucurrentgroup[1])
    scatter3 = ax3.scatter(dftest_embedding[0], dftest_embedding[1],s=150,marker='*')
    fig3=plt.figure(figsize=(50,95))

    for i in range(len(dfcurrentgroup)):
        ax3.annotate(dfcurrentgroup['DNA sample ethnicity'][i], (dfucurrentgroup[0][i], dfucurrentgroup[1][i]))
    st.pyplot()
