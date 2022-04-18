# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 17:59:38 2020

@author: Y.H. Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ discription ==============================####
#
#================================== input =====================================
#
#================================== output ====================================
#
#================================ parameters ==================================
#
#================================== example =================================== 
#
#================================== warning ===================================
# seems the length of head can only be 3
####=======================================================================####
"""
import pandas as pd

def abund_table2lefse(df,head_index,head_data):
    head_dict=dict(zip(head_index,head_data))
    head_df=pd.DataFrame.from_dict(head_dict,orient='index')
    
    lefse_df=pd.concat([head_df,df],axis=0)
    
    return lefse_df

def rel_abun2lefse(dp,Class):
    basename=os.path.basename(dp)
    
    df=pd.read_csv(dp,sep=",",index_col=0).T
    
    subClass=['subClass']*len(Class)
    #Subject=["305A_07022019","305A_07302019","396_07022019","396_07302019","473B_07022019","473B_07302019","705E_07012019","705E_07302019","7108E_07022019","7108E_07302019","711E_07012019","711E_07302019","778E_07012019","778E_07302019","782E_07012019","782E_07302019"]
    Subject=list(df.columns)
    df.columns=range(df.shape[1])
    
    """
    Class=['DH','DH','DH','HM','HM','HM','DH','DH',\
              'DH','DH','DH','DH','HM','HM','HM','HM','HM','HM']
    
    subClass=['hinny','hinny','hinny','mule','mule','mule',\
                 'donkey','donkey','donkey','donkey',\
                 'donkey','donkey','horse','horse','horse','horse','horse','horse']
    Subject=["hinny3110","hinny2611","hinny3742","mule3163","mule31631","mule4285","donkey3611","donkey3446",\
        "donkey3418","donkey4058","donkey4282","donkey4687","horse3474","horse4964","horse3912","horse3729",\
        "horse4231","horse4986"]
    """
    """
    Class=['donkey+hinny','donkey+hinny','donkey+hinny','horse+mule','horse+mule','horse+mule','donkey+hinny','donkey+hinny','donkey+hinny',\
           'donkey+hinny','donkey+hinny','donkey+hinny','horse+mule','horse+mule','horse+mule','horse+mule','horse+mule','horse+mule']
    subClass=['hmdh']*18
    Subject=["hinny_3110","hinny_2611","hinny_3742","mule_3163","mule_31632","mule_4285","donkey_3611","donkey_3446","donkey_3418","donkey_4058","donkey_4282","donkey_4687","horse_3474","horse_4964","horse_3912","horse_3729","horse_4231","horse_4986"]
    """
    head_index=['Class','subClass','Subject']
    head_data=[Class,subClass,Subject]
    
    lefse_df=abund_table2lefse(df,head_index,head_data)
    
    lefse_df.to_csv("../lefse/{}.lefse.tsv".format(basename),header=None,sep='\t')
        


if __name__=='__main__':
    import os
    if not os.path.exists('../lefse/'):
        os.makedirs('../lefse/')
    """
    for i in range(1,9):
        dp=r"../taxa_abun/rel_abun/sample12_rel_abun.{}.rmU.euk.csv".format(str(i))
        Class=['donkey','donkey','donkey','donkey','donkey','donkey',\
               'horse','horse','horse','horse','horse','horse']
        rel_abun2lefse(dp,Class)
    """
    for i in ['ko','eggnog','level4ec','pfam','rxn']:
        #dp=r'../metabolism/humann3/horsedonkey_genefamilies_uniref90names_cpm_{}_unstratified.named.tsv.rel_abun.csv'.format(i)
        dp=r'C:\CurrentProjects\equids_MHC\Prj1\humann3\horsemuledonkeyhinny_genefamilies_uniref90names_cpm_{}_unstratified.named.tsv.rela_abun_style.csv'.format(i)
        """
        Class=['donkey','donkey','donkey','donkey','donkey','donkey',\
               'horse','horse','horse','horse','horse','horse']
        """
        Class=['donkey+hinny']*9+['horse+mule']*9
        rel_abun2lefse(dp,Class)
    
    
