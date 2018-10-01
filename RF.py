import numpy as np
from sklearn.metrics import accuracy_score
import copy
import random

def corrona(drug):

    import pandas as pd

    path = '/Users/asherameli/Desktop/RA_response_masterfile.csv'
    df = pd.read_csv(path, header=None,
                       names=[
        'SUBJECT','Row.names','Sample.ID','Suffix','Specimen.Name','Instrument','PDF','X28S.18S',
      'RIN','Qubit.ug.ul','Total..ug.','DETD','STUDY.ID','SHIPMENT.NUMBER',	'SAMPLE.ID',
      'CONTAINER',	'ACCESSION.NUMBER',	'DRAW.DATE','DRAW.TIME','SITE...',	'INITIALS',
      'DATE.OF.BIRTH','GENDER',	'VISIT','QC',
    'Scipher_id',	'grp','visit','visitdate','CDate',
      'age','gender','duration_ra','race_grp','race_hispanic','hx_fib','hx_CVD','hx_ser_inf',
      'hx_anycancer','bmi',	'newsmoker','dose_freq_bio','con_pred',	'dose_pred','con_mtx',
      'dose_mtx', 'con_cdmard',	'dose_cdmard','usresultsCRP','das28crp','swollen_jts_28','tender_jts_28',
      'pt_global_assess','md_global_assess','cdai','di','pt_pain','eular_dascrp','status_6m',
      'reason_disc_bio','switch','ccpposever_new','usresultsCCP3','rfposever_new','usresultsRF',
        'eular_das28crp','mcid', 'endpoint_cdai','sdai','endpoint_sdai',
          'endpoint_das28crp', 'acr_20', 'acr_50', 'acr_70'
                           ,'vote','decision_vote'
                       ])


    user2feature_base = {}
    user2feature3 = {}
    user2feature6 = {}

    for index , data in df.iterrows():

        if index==0: continue
        user_id =data['Scipher_id']

        if drug =='all':

            if data['visit'] == '0': user2feature_base[user_id] = data
            if data['visit'] == '3': user2feature3[user_id] = data
            if data['visit'] == '6': user2feature6[user_id] = data


        if drug == 'enbrel':
            if data['grp'] ==drug:

                if data['visit'] == '0': user2feature_base[user_id] = data
                if data['visit'] == '3': user2feature3[user_id] = data
                if data['visit'] == '6': user2feature6[user_id] = data

    return [user2feature_base,user2feature3,user2feature6]


features=[
    'age','duration_ra','hx_fib','hx_CVD','hx_ser_inf',
      'hx_anycancer','bmi',
      'dose_mtx', 'con_cdmard',	'dose_cdmard','usresultsCRP','das28crp','swollen_jts_28','tender_jts_28',
      'pt_global_assess','md_global_assess','cdai','di','pt_pain'
]

all_responses = [ 'mcid', 'endpoint_cdai','endpoint_das28crp','eular_das28crp', 'acr_20', 'acr_50', 'acr_70']


from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split


condition=['3','all']
title = condition[0]+' months; '+condition[1]+' drugs'


corrona_data = corrona(condition[1])
user2feature_base = corrona_data[0]
user2feature_3 = corrona_data[1]
user2feature_6 = corrona_data[2]

if condition[0]=='3': user2feature = user2feature_3
if condition[0]=='6': user2feature = user2feature_6


gender_2_index = {'female':0,'male':1}
race_grp_2_index={'Mixed Race':0,'Other':1,'Black':2,'Native American':3,'White':4,'Asian':5}
race_hispanic_2_index={'not':0,'hispanic':1}
newsmoker_2_index={'previous':0,'current':1,'never':2}
con_pred_2_index={'no':0,'yes':1}
response_2_index = {'poor':0,'moderate':0,'good':1}
response_2_index2 = {'non-responder':0,'responder':1}

accuracy=[]
accuracy_rnd=[]
for response_choice in all_responses:

    X = []
    y =[]
    for user, data in user2feature.iteritems():

        data_base = user2feature_base[user]
        info = [data_base[f] for f in features]

        info_clean=[]
        for f in info:

            if f in gender_2_index.keys(): f = gender_2_index[f]
            if f in race_grp_2_index.keys(): f= race_grp_2_index[f]
            if f in race_hispanic_2_index.keys(): f=  race_hispanic_2_index[f]
            if f in newsmoker_2_index.keys(): f= newsmoker_2_index[f]
            if f in con_pred_2_index.keys(): f= con_pred_2_index[f]

            info_clean.append(float(f))

        response = data[response_choice]

        if response in response_2_index.keys(): response = response_2_index[response]
        if response in response_2_index2.keys(): response = response_2_index2[response]

        if np.any(np.isnan(info_clean)) ==True: continue
        if np.any(np.isnan(response)) == True: continue

        X.append(info_clean)
        y.append(float((response)))


    print response_choice, len(y),np.sum(y), np.sum(y)/len(y)

    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state = 0)
    forest = RandomForestClassifier(n_estimators=10, random_state=0)
    forest.fit(X_train,y_train)


    y_pred = forest.predict(X_test)
    acc =  accuracy_score(y_test, y_pred)
    accuracy.append(acc)
    #-------------------------------
    yy = copy.copy(y)

    random_acc=[]
    for i in xrange(1000):
        yy= random.sample(yy, len(yy))
        X_train, X_test, y_train, y_test = train_test_split(X, yy, random_state=0)

        forest = RandomForestClassifier(n_estimators=10, random_state=0)
        forest.fit(X_train, y_train)
        y_pred = forest.predict(X_test)
        acc_rnd = accuracy_score(y_test, y_pred)
        random_acc.append(acc_rnd)
    accuracy_rnd.append(acc_rnd)
    # -------------------------------
    print np.mean(random_acc),acc
    print '-------------------'


from matplotlib import pyplot as plt
plt.style.use('ggplot')
plt.rcParams['text.usetex'] = True

plt.scatter(xrange(len(accuracy)),accuracy, s=40, c='r', label='actual')
plt.scatter(xrange(len(accuracy)),accuracy_rnd, s=40,c='b', label='random')
plt.legend(loc='lower right')
plt.title(title)
plt.ylim([0,1])
plt.ylabel('Accuracy')
plt.xticks(xrange(len(accuracy)),[res.replace('_',' - ').upper() for res in all_responses], rotation=80, fontsize=12)
font = {'family': 'Helvetica', 'weight': 'bold', 'size': 20}
plt.rc('font', **font)
plt.savefig('out/ML/RF/'+title+'.png', dpi=200, bbox_inches='tight'); plt.close()





