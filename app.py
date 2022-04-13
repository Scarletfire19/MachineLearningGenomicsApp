import streamlit as st
from threading import activeCount
import matplotlib.pyplot as plt
import pickle
import umap
import io
import numpy as np
import pandas as pd
import streamlit_authenticator as stauth
import os
import plotly.express as px



st.set_option('deprecation.showPyplotGlobalUse', False)

#https://wallpapercave.com/wp/wp4594030.jpg
#https://wallpapercave.com/wp/wp2170040.jpg
#https://blog.23andme.com/wp-content/uploads/2008/09/novembreblogpostfig.jpg
#https://media.istockphoto.com/photos/background-picture-id185314478?b=1&k=20&m=185314478&s=612x612&w=0&h=XTWvVxe3IdiV6B5bJDi2DAICfqmyaZxWnWyXY_Lw5kk=

st.title("Machine Learning Genomics App")

names = ['Blue Bear','Cream Rabbit']
usernames = ['bhalu','khwabekhargosh']
passwords = ['bear','rabbit']
hashed_passwords = stauth.Hasher(passwords).generate()
authenticator = stauth.Authenticate(names,usernames,hashed_passwords,'some_cookie_name','some_signature_key',cookie_expiry_days=30)
name, authentication_status, username = authenticator.login('Login','main')

if authentication_status:
    st.write('Welcome *%s*' % (name))
    st.title('Some content')
elif authentication_status == False:
    st.error('Username/password is incorrect')
elif authentication_status == None:
    st.warning('Please enter your username and password')

if authentication_status:
    authenticator.logout('Logout', 'main')
    st.write('Welcome *%s*' % (name))
    st.title('Some content')
elif authentication_status == False:
    st.error('Username/password is incorrect')


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

dfumap=pd.read_csv('umap.csv')
dftsne=pd.read_csv('tSNE.csv')


uploaded_file = st.file_uploader("Enter G25 co-ordinates")
if uploaded_file is not None:
     input = pd.read_csv(uploaded_file)
     st.write(input)
      
          
inputhaplogroup = st.text_input("Enter Y Haplogroup")

def euclidean_distance(w, q):
    n = 25 
    return sum([(w[i] - q[i]) ** 2 for i in range(n)]) ** 0.5

p=Xcombined.iloc[735:740]
dfdistances=dfcombined
distances=[]



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
k=0


#dfdistances
#print(dfdistances['DNA sample ethnicity and id'].iloc[:3])

Tools = st.selectbox("Choose your Tool", ["Genetic Distance Tool", "PCA(Principal Component Analysis) Tool","ML Ancestry Tool","Ancient DNA Lineage Tool","tSNE","Umap"]) 

if Tools == "Genetic Distance Tool":
     page_bg_img = '''
     <style>
     body {
     background-image: url("https://blog.23andme.com/wp-content/uploads/2008/09/novembreblogpostfig.jpg");
     background-size: cover;
     }
     </style>
     '''
     st.markdown(page_bg_img, unsafe_allow_html=True)


     st.title("Genetic Distance Tool")
     Xinput=input.drop(columns=['DNA sample ethnicity and id'])


     for i in range(len(Xcombined)):
          distances.append(euclidean_distance(Xcombined.iloc[i],Xinput.iloc[3]))
     dfdistances['distances']=distances

     dfdistances=dfdistances.sort_values(by=['distances'])
     st.dataframe(dfdistances)


elif Tools == "ML Ancestry Tool":
    
     page_bg_img2 = '''
     <style>
     body {
     background-image: url("https://cdn-prod.medicalnewstoday.com/content/images/articles/318/318492/dna-strands.jpg");
     background-size: cover;
     }
     </style>
     '''
     st.markdown(page_bg_img2, unsafe_allow_html=True)
     

     st.title("ML Ancestry Tool")
     st.dataframe(ancestry)


elif Tools == "PCA(Principal Component Analysis) Tool":
    
     st.title("PCA(Principal Component Analysis) Tool")
        
     pd.options.plotting.backend = "plotly"
     figplotly0 = dfcurrentgroup.plot.scatter(x="1", y="2",text="DNA sample ethnicity")
     figplotly0 = input.plot.scatter(x="1", y="2",text="DNA sample ethnicity")
     st.plotly_chart(figplotly0)   

     #fig, ax = plt.subplots(figsize=(45, 99))
     #ax.scatter(dfcurrentgroup['1'], dfcurrentgroup['2'],s = 1)
     #ax.scatter(input['1'], input['2'],s = 5)
     #ax.scatter(point['1'],point['2'],s=500)

     #for i in range(len(dfcurrentgroup)):
     #     ax.annotate(dfcurrentgroup['DNA sample ethnicity'][i], (dfcurrentgroup['1'][i], dfcurrentgroup['2'][i]))
     #st.pyplot()


elif Tools == "Ancient DNA Lineage Tool":
     page_bg_img3 = '''
     <style>
     body {
     background-image: url("https://blog.23andme.com/wp-content/uploads/2008/09/novembreblogpostfig.jpg");
     background-size: cover;
     }
     </style>
     '''
     st.markdown(page_bg_img3, unsafe_allow_html=True)

    
     ancienthaplogroup = dfancienthpg[dfancienthpg['Assigned Mutation'].str.contains(inputhaplogroup)]     
     st.title("Ancient DNA Lineage Tool")

     fig2, ax2 = plt.subplots(figsize=(30, 24))
     ax2.scatter(dfancienthpg['1'], dfancienthpg['2'],s = 1)
     ax2.scatter(input['1'], input['2'],s = 5)
     ax2.scatter(ancienthaplogroup['1'], ancienthaplogroup['2'],s = 100)

     #ax.scatter(point['1'],point['2'],s=500)
     for i in range(len(dfancienthpg)):
          ax2.annotate(dfancienthpg['Assigned Mutation'][i], (dfancienthpg['1'][i], dfancienthpg['2'][i]))
     st.pyplot()
    
elif Tools == "tSNE":
     pd.options.plotting.backend = "plotly"
     figplotly1 = dftsne.plot.scatter(x="0", y="1",text="DNA sample ethnicity")
     st.plotly_chart(figplotly1, use_container_width=True)
    
elif Tools == "Umap":
     pd.options.plotting.backend = "plotly"
     figplotly2 = dfumap.plot.scatter(x="0", y="1",text="DNA sample ethnicity")
     st.plotly_chart(figplotly2, use_container_width=True)
    


import os
#filename = 'umap_model.sav'
#current_path=os.getcwd()
#modelumap_path = os.path.join(current_path, 'umap_model.sav')
#loaded_model = pickle.load(open(modelumap_path, 'rb'))

pickle_in = open('umap_model.pkl', 'rb') 
loaded_model = pickle.load(pickle_in)

test_embedding = loaded_model.transform(input)
dftest_embedding=pd.DataFrame(test_embedding)

fig5, ax5 = plt.subplots(figsize=(50, 200))
scatter = ax5.scatter(dftest_embedding[0], dftest_embedding[1],s=120,marker='*')
fig=plt.figure(figsize=(55,90))
st.pyplot()

