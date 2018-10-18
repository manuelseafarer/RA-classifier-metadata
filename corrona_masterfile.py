import pandas as pd

def eular_das28crp(data_base,data):

    das28crp_endpoint = float(data['das28crp'])
    das28crp_baseline = float(data_base['das28crp'])
    improvement = das28crp_baseline - das28crp_endpoint

    if das28crp_endpoint <= 2.7:

        if improvement > 1.2: response='good'
        if 0.6 <improvement <= 1.2: response = 'moderate'
        if improvement <= 0.6: response = 'poor'

    if 2.7 <das28crp_endpoint <= 4.1:

        if improvement > 1.2: response = 'moderate'
        if 0.6 < improvement <= 1.2: response = 'moderate'
        if improvement <= 0.6: response = 'poor'

    if das28crp_endpoint > 4.1:

        if improvement > 1.2: response = 'moderate'
        if 0.6 < improvement <= 1.2: response = 'poor'
        if improvement <= 0.6: response = 'poor'

    return response

def mcid(data_base, data):

    cdai = float(data['cdai'])
    cdai_base = float(data_base['cdai'])
    improvement = cdai_base - cdai

    response = 'non-responder'
    if 10.1<= cdai_base<22 and improvement>6: response='responder'
    if cdai_base >22 and improvement>12: response = 'responder'

    return response

def mcid20(data_base, data):

    cdai = float(data['cdai'])
    cdai_base = float(data_base['cdai'])
    improvement = cdai_base - cdai

    response = 'non-responder'
    if 10.1<= cdai_base<20 and improvement>6: response='responder'
    if cdai_base >20 and improvement>12: response = 'responder'

    return response

def endpoint_cdai(data):

    cdai = float(data['cdai'])
    response ='non-responder'
    if 0<=cdai<=10.0: response='responder'
    return response

def endpoint_sdai(data):

    #TJC28 + SJC28 + patient global score + physician global score + CRP(mg/dl)
    tender_jts_28 = float(data['tender_jts_28'])
    swollen_jts_28 = float(data['swollen_jts_28'])
    pt_global_assess = float(data['pt_global_assess'])
    md_global_assess = float(data['md_global_assess'])
    usresultsCRP = float(data['usresultsCRP'])
    sdai = tender_jts_28 + swollen_jts_28 + pt_global_assess + md_global_assess + usresultsCRP
    response ='non-responder'
    if 0<=sdai<=11.0: response='responder'
    return [sdai, response]


def endpoint_das28crp(data):

    das28crp = float(data['das28crp'])
    response ='non-responder'
    if 0<=das28crp<2.9: response='responder'

    return response

def acr(data_base,data,intensity):

    main_switch='off'
    tender_jts_28_base = float(data_base['tender_jts_28'])
    tender_jts_28 = float(data['tender_jts_28'])
    try:
        main_improvement1 = (tender_jts_28_base - tender_jts_28)/tender_jts_28_base
    except ZeroDivisionError:
        main_improvement1 = 0.0

    swollen_jts_28_base = float(data_base['swollen_jts_28'])
    swollen_jts_28 = float(data['swollen_jts_28'])
    try:
        main_improvement2 = (swollen_jts_28_base - swollen_jts_28) / swollen_jts_28_base
    except ZeroDivisionError:
        main_improvement2 = 0.0



    if main_improvement1 >= intensity or main_improvement2>= intensity: main_switch='on'

    secondary_switch = 0
    di_base = float(data_base['di'])
    di = float(data['di'])
    try:
        di_improvement = (di_base - di) / di_base
    except ZeroDivisionError:
        di_improvement = 0.0
    if di_improvement> intensity: secondary_switch=secondary_switch+1

    pt_pain_base = float(data_base['pt_pain'])
    pt_pain = float(data['pt_pain'])
    try:
        pt_pain_improvement = (pt_pain_base - pt_pain) / pt_pain_base
    except ZeroDivisionError:
        pt_pain_improvement = 0.0
    if pt_pain_improvement > intensity: secondary_switch = secondary_switch + 1


    pt_global_assess_base = float(data_base['pt_global_assess'])
    pt_global_assess = float(data['pt_global_assess'])
    try:
        pt_global_assess_improvement = (pt_global_assess_base - pt_global_assess) / pt_global_assess_base
    except ZeroDivisionError:
        pt_global_assess_improvement=0.0
    if pt_global_assess_improvement > intensity: secondary_switch = secondary_switch + 1

    md_global_assess_base = float(data_base['md_global_assess'])
    md_global_assess = float(data['md_global_assess'])
    md_global_assess_improvement = (md_global_assess_base - md_global_assess) / md_global_assess_base
    if md_global_assess_improvement > intensity: secondary_switch = secondary_switch + 1

    usresultsCRP_base = float(data_base['usresultsCRP'])
    usresultsCRP = float(data['usresultsCRP'])
    usresultsCRP_improvement = (usresultsCRP_base - usresultsCRP) / usresultsCRP_base
    if usresultsCRP_improvement > intensity: secondary_switch = secondary_switch + 1

    response= 'non-responder'
    if main_switch =='on' and secondary_switch>=3:response='responder'

    return response


original_header=['SUBJECT','Row.names','Sample.ID','Suffix','Specimen.Name','Instrument','PDF','X28S.18S',
                              'RIN','Qubit.ug.ul','Total..ug.','DETD','STUDY.ID','SHIPMENT.NUMBER',	'SAMPLE.ID',
                              'CONTAINER',	'ACCESSION.NUMBER',	'DRAW.DATE','DRAW.TIME','SITE...',	'INITIALS',
                              'DATE.OF.BIRTH','GENDER',	'VISIT','QC','Scipher_id',	'grp','visit','visitdate','CDate',
                              'age','gender','duration_ra','race_grp','race_hispanic','hx_fib','hx_CVD','hx_ser_inf',
                              'hx_anycancer','bmi',	'newsmoker','dose_freq_bio','con_pred',	'dose_pred','con_mtx',
                              'dose_mtx', 'con_cdmard',	'dose_cdmard','usresultsCRP','das28crp','swollen_jts_28','tender_jts_28',
                              'pt_global_assess','md_global_assess','cdai','di','pt_pain','eular_dascrp','status_6m',
                              'reason_disc_bio','switch','ccpposever_new','usresultsCCP3','rfposever_new','usresultsRF']
def corrona():

    path = 'data/RA_ClinialCovariates.csv'
    df = pd.read_csv(path, header=None,
                       names=original_header)

    user2baseline = {}
    user2visit3 = {}
    user2visit6 = {}

    for index , data in df.iterrows():

        scipher_id =data['Scipher_id']

        if data['visit']=='0': user2baseline[scipher_id] = data
        if data['visit']=='3': user2visit3[scipher_id] = data
        if data['visit']=='6': user2visit6[scipher_id] = data

    return [user2baseline,user2visit3,user2visit6]

corrona_data = corrona()
user2baseline = corrona_data[0]
user2visit3 = corrona_data[1]
user2visit6 = corrona_data[2]
all_users3=user2visit3.keys()
all_users6=user2visit6.keys()



eular_das28crp_dic3={}
mcid_dic3={}
mcid20_dic3={}
acr_dic3_20={}
acr_dic3_50={}
acr_dic3_70={}
endpoint_cdai_dic3={}
endpoint_das_dic3={}
endpoint_sdai_dic3={}

for user in all_users3:

    data = user2visit3[user]
    data_base = user2baseline[user]
    response_eular = eular_das28crp(data_base,data)
    eular_das28crp_dic3[user] = response_eular

    response_mcid = mcid(data_base, data)
    mcid_dic3[user] = response_mcid

    response_mcid20 = mcid20(data_base, data)
    mcid20_dic3[user] = response_mcid20

    response_acr20 = acr(data_base, data, 0.2)
    response_acr50 = acr(data_base, data, 0.5)
    response_acr70 = acr(data_base, data, 0.7)

    acr_dic3_20[user]=response_acr20
    acr_dic3_50[user] = response_acr50
    acr_dic3_70[user] = response_acr70


    response_end_cdi=endpoint_cdai(data)
    response_end_das = endpoint_das28crp(data)
    response_end_sdai = endpoint_sdai(data)[1]

    endpoint_cdai_dic3[user] = response_end_cdi
    endpoint_das_dic3[user] = response_end_das

    endpoint_sdai_dic3[user] = response_end_sdai

eular_das28crp_dic6 = {}
mcid_dic6 = {}
mcid20_dic6 = {}
acr_dic6_20 = {}
acr_dic6_50 = {}
acr_dic6_70 = {}
endpoint_cdai_dic6 = {}
endpoint_das_dic6 = {}
endpoint_sdai_dic6 = {}

for user in all_users6:
    data = user2visit6[user]
    data_base = user2baseline[user]
    try:
        response_eular = eular_das28crp(data_base, data)

    except UnboundLocalError:
        response_eular='NA'
    eular_das28crp_dic6[user] = response_eular

    response_mcid = mcid(data_base, data)
    mcid_dic6[user] = response_mcid

    response_mcid20 = mcid20(data_base, data)
    mcid20_dic6[user] = response_mcid20

    response_acr20 = acr(data_base, data, 0.2)
    response_acr50 = acr(data_base, data, 0.5)
    response_acr70 = acr(data_base, data, 0.7)

    acr_dic6_20[user] = response_acr20
    acr_dic6_50[user] = response_acr50
    acr_dic6_70[user] = response_acr70

    response_end_cdi = endpoint_cdai(data)
    response_end_das = endpoint_das28crp(data)
    response_end_sdai = endpoint_sdai(data)[1]

    endpoint_cdai_dic6[user] = response_end_cdi
    endpoint_das_dic6[user] = response_end_das
    endpoint_sdai_dic6[user] = response_end_sdai



print_header = ('SUBJECT','Row.names','Sample.ID','Suffix','Specimen.Name','Instrument','PDF','X28S.18S',
      'RIN','Qubit.ug.ul','Total..ug.','DETD','STUDY.ID','SHIPMENT.NUMBER',	'SAMPLE.ID',
      'CONTAINER',	'ACCESSION.NUMBER',	'DRAW.DATE','DRAW.TIME','SITE...',	'INITIALS',
      'DATE.OF.BIRTH','GENDER',	'VISIT','QC','Scipher_id',	'grp','visit','visitdate','CDate',
      'age','duration_ra','race_grp','race_hispanic','hx_fib','hx_CVD','hx_ser_inf',
      'hx_anycancer','bmi',	'newsmoker','dose_freq_bio','con_pred',	'dose_pred','con_mtx',
      'dose_mtx', 'con_cdmard',	'dose_cdmard','usresultsCRP','das28crp','swollen_jts_28','tender_jts_28',
      'pt_global_assess','md_global_assess','cdai','di','pt_pain','eular_dascrp','status_6m',
      'reason_disc_bio','switch','ccpposever_new','usresultsCCP3','rfposever_new','usresultsRF',
        'eular_das28crp', 'mcid','mcid_20', 'endpoint_cdai','endpoint_das28crp', 'acr_20', 'acr_50', 'acr_70','acr50-cdai-10')



f = open('data/RA_response_masterfile_additional.csv', 'w')
f.write(','.join('{}'.format(val) for val in print_header)); f.write('\n')



path = 'data/RA_ClinialCovariates.csv'
df = pd.read_csv(path, header=None,names=original_header)



positive_classes = ['responder', 'good']
for index , data in df.iterrows():

    if index ==0: continue

    scipher_id = data['Scipher_id']

    sdai = endpoint_sdai(data)[0]

    if data['visit'] == '0':

        eular_das28crp_score = ' '
        mcid_score = ' '
        mcid20_score = ' '
        endpoint_cdai_score = ' '
        endpoint_sdai_score = ' '
        endpoint_das28crp_score = ' '
        acr_20_score = ' '
        acr_50_score = ' '
        acr_70_score = ' '
        acr50_cdai_10_score= ' '

        '''
        decision= ' '
        counter = ' '
        '''

    if data['visit'] == '3':
        counter = 0
        eular_das28crp_score = eular_das28crp_dic3[scipher_id]
        if eular_das28crp_score in positive_classes: counter= counter+1

        mcid_score = mcid_dic3[scipher_id]
        if mcid_score in positive_classes: counter= counter+1

        mcid20_score = mcid20_dic3[scipher_id]
        if mcid20_score in positive_classes: counter = counter + 1

        endpoint_cdai_score = endpoint_cdai_dic3[scipher_id]
        if endpoint_cdai_score in positive_classes: counter= counter+1

        endpoint_sdai_score = endpoint_sdai_dic3[scipher_id]
        if endpoint_sdai_score in positive_classes: counter = counter+1

        endpoint_das28crp_score = endpoint_das_dic3[scipher_id]
        if endpoint_das28crp_score in positive_classes: counter =  counter+1

        acr_20_score = acr_dic3_20[scipher_id]
        if acr_20_score in positive_classes: counter = counter+1

        acr_50_score = acr_dic3_50[scipher_id]

        if acr_50_score in positive_classes: counter = counter+1

        acr_70_score = acr_dic3_70[scipher_id]
        if acr_70_score in positive_classes: counter = counter+1

        acr50_cdai_10_score='non-responder'
        if (acr_50_score in positive_classes) or (endpoint_cdai_score in positive_classes): acr50_cdai_10_score='responder'

        '''
        if counter>4: decision ='responder'
        if counter<=4: decision = 'non-responder'
        '''

    if data['visit'] == '6':

        counter = 0
        eular_das28crp_score = eular_das28crp_dic6[scipher_id]
        if eular_das28crp_score in positive_classes: counter = counter + 1

        mcid_score = mcid_dic6[scipher_id]
        if mcid_score in positive_classes: counter = counter + 1

        mcid20_score = mcid20_dic6[scipher_id]
        if mcid20_score in positive_classes: counter = counter + 1

        endpoint_cdai_score = endpoint_cdai_dic6[scipher_id]
        if endpoint_cdai_score in positive_classes: counter = counter + 1

        endpoint_sdai_score = endpoint_sdai_dic6[scipher_id]
        if endpoint_sdai_score in positive_classes: counter = counter + 1

        endpoint_das28crp_score = endpoint_das_dic6[scipher_id]
        if endpoint_das28crp_score in positive_classes: counter = counter + 1

        acr_20_score = acr_dic6_20[scipher_id]
        if acr_20_score in positive_classes: counter = counter + 1

        acr_50_score = acr_dic6_50[scipher_id]
        if acr_50_score in positive_classes: counter = counter + 1

        acr_70_score = acr_dic6_70[scipher_id]
        if acr_70_score in positive_classes: counter = counter + 1

        acr50_cdai_10_score = 'non-responder'
        if (acr_50_score in positive_classes) or (endpoint_cdai_score in positive_classes): acr50_cdai_10_score = 'responder'
        '''
        if counter > 4: decision = 'responder'
        if counter <= 4: decision = 'non-responder'
        '''



    print_data = [data['SUBJECT'],data['Row.names'],data['Sample.ID'],data['Suffix'],data['Specimen.Name'],data['Instrument'],
             data['PDF'],data['X28S.18S'],data['RIN'],data['Qubit.ug.ul'],data['Total..ug.'],data['DETD'],data['STUDY.ID'],
             data['SHIPMENT.NUMBER'],data['SAMPLE.ID'],data['CONTAINER'],data['ACCESSION.NUMBER'],data['DRAW.DATE'],'DRAW.TIME','SITE...',	'INITIALS',
            data['DATE.OF.BIRTH'],data['GENDER'],data['VISIT'],data['QC'],data['Scipher_id'],data['grp'],data['visit'],
             data['visitdate'],data['CDate'],data['age'],data['duration_ra'],data['race_grp'],data['race_hispanic'],
             data['hx_fib'],data['hx_CVD'],data['hx_ser_inf'],data['hx_anycancer'],data['bmi'],	data['newsmoker'],
             data['dose_freq_bio'],data['con_pred'],data['dose_pred'],data['con_mtx'],data['dose_mtx'], data['con_cdmard'],
             data['dose_cdmard'],data['usresultsCRP'],data['das28crp'],data['swollen_jts_28'],data['tender_jts_28'],
            data['pt_global_assess'],data['md_global_assess'],data['cdai'],data['di'],data['pt_pain'],data['eular_dascrp'],
             data['status_6m'],data['reason_disc_bio'],data['switch'],data['ccpposever_new'],data['usresultsCCP3'],
             data['rfposever_new'],data['usresultsRF'],eular_das28crp_score,mcid_score,mcid20_score, endpoint_cdai_score,
             endpoint_das28crp_score,acr_20_score,acr_50_score,acr_70_score, acr50_cdai_10_score]

    f.write(','.join('{}'.format(val) for val in print_data))
    f.write('\n')











