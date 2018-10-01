from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import random
import copy
import numpy as np
from matplotlib import pyplot as plt
plt.style.use('ggplot')
plt.rcParams['text.usetex'] = True

#-----------
import numpy as np

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
#response_choice= all_responses[6]



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
#----------
for response_choice in all_responses:

    X = []
    y =[]

    for user, data in user2feature.iteritems():

        data_base = user2feature_base[user]
        info = [data_base[f] for f in features]

        info_clean = []
        for f in info:

            if f in gender_2_index.keys(): f = gender_2_index[f]
            if f in race_grp_2_index.keys(): f = race_grp_2_index[f]
            if f in race_hispanic_2_index.keys(): f = race_hispanic_2_index[f]
            if f in newsmoker_2_index.keys(): f = newsmoker_2_index[f]
            if f in con_pred_2_index.keys(): f = con_pred_2_index[f]

            info_clean.append(float(f))

        response = data[response_choice]

        if response in response_2_index.keys(): response = response_2_index[response]
        if response in response_2_index2.keys(): response = response_2_index2[response]

        if np.any(np.isnan(info_clean)) == True: continue
        if np.any(np.isnan(response)) == True: continue

        X.append(info_clean)
        y.append(float((response)))



    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state = 0)
    forest = RandomForestClassifier(n_estimators=10, random_state=0)
    forest.fit(X_train,y_train)

    y_pred = forest.predict(X_test)
    accuracy = float(accuracy_score(y_test, y_pred))
    print accuracy

    #---------------------------




    feature2_impact={}
    for id_f, f in enumerate(features):

        target_vals = []
        for r in X:

            target_vals.append(r[id_f])

        temp=[]
        XX = copy.copy(X)


        for i in xrange(1000):

            random_target_vals = random.sample(target_vals, len(target_vals))
            for id, val in enumerate(XX):
                XX[id][id_f] = random_target_vals[id]

            X_train, X_test, y_train, y_test = train_test_split(XX, y, random_state=0)
            forest = RandomForestClassifier(n_estimators=10, random_state=0)
            forest.fit(X_train, y_train)

            y_pred = forest.predict(X_test)



            acc = float(accuracy_score(y_test, y_pred))

            temp.append(float(acc))
        feature2_impact[f]= (accuracy - np.mean(temp))/np.std(temp)


    features_importance=[]
    for f in features:
        features_importance.append(feature2_impact[f])


    print features_importance
    title = condition[0] + ' months; ' + condition[1] + ' drugs; ' + response_choice.replace('_', '- ').upper()

    plt.xlim([-0.5,len(features)])
    plt.bar(range(len(features_importance)), features_importance, align='center')
    plt.xticks(xrange(len(features)),[res.replace('_',' - ').upper() for res in features], rotation=80, fontsize=10)
    plt.title(title, fontsize=15)
    plt.ylabel('Z-score')
    font = {'family': 'Helvetica', 'weight': 'bold', 'size': 16}
    plt.rc('font', **font)
    plt.savefig('out/ML/RF/'+title+' - '+response_choice+'.png', dpi=200, bbox_inches='tight'); plt.close()