# -*- coding: utf-8 -*-
"""
Created on Wed May 13 11:24:01 2020

@author: Y.H. Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ discription ==============================####
## part 1: nonparametric test and output
## part 2: give the genes list which significant changed
## part 2 need to use interparameter in part 1
#================================== input =====================================
#
#================================== output ====================================
#
#================================ parameters ==================================
#
#================================== example ===================================

#================================== warning ===================================
#
####=======================================================================####
"""

from scipy import stats
import pandas as pd
import sys

def extract_lists(df,column,index1,index2):
    list1=df.loc[index1,column].tolist()
    list2=df.loc[index2,column].tolist()
    return [list1,list2]

def normalty_Wilcoxon_dual_group_test(df,group_1,group_2,paired=False,only_two_sided=True):
    # lists of abundance, list1 and list2 responds to group1 and group2
    lists=[] #[[list1,list2],[]....]
    normalty=[]
    normalty_p=[]
    Wilcoxon=[]
    Wilcoxon_p=[]
    taxa_list=[] #["Firm","Bact"...]
    for col in df.columns:
        if sum(df[col])==0:
            print(col," is zero for all samples. skip.")
            continue
        taxa_list.append(col)
        [list1,list2]=extract_lists(df,col,group_1,group_2)
        list12=list1+list2
        lists.append([list1,list2])

        alpha=1e-3
        k, p = stats.kstest(list12,'norm')
        if p<alpha:
            normalty.append("N")
            normalty_p.append(p)
        else:
            normalty.append("Y")
            normalty_p.append(p)

        alpha=0.05
        if only_two_sided==False:
            if paired==True:
                k,p=stats.wilcoxon(list1,list2,alternative="less")
            else:
                k,p=stats.mannwhitneyu(list1,list2,alternative="less")
            if p<alpha:
                Wilcoxon.append("increase")
                Wilcoxon_p.append(p)
            else:
                # print("less_p: ",p)
                if paired==True:
                    k,p=stats.wilcoxon(list1,list2,alternative="greater")
                else:
                    k,p=stats.mannwhitneyu(list1,list2,alternative="greater")
                if p<alpha:
                    Wilcoxon.append("decrease")
                    Wilcoxon_p.append(p)
                else:
                    if paired==True:
                        k,p=stats.wilcoxon(list1,list2,alternative="two-sided")
                    else:
                        k,p=stats.mannwhitneyu(list1,list2,alternative="two-sided")
                    Wilcoxon.append("equal")
                    Wilcoxon_p.append(p)
        else:
            if paired==True:
                k,p=stats.wilcoxon(list1,list2,alternative="two-sided")
            else:
                k,p=stats.mannwhitneyu(list1,list2,alternative="two-sided")
            if p<alpha:
                Wilcoxon.append("unequal")
            else:
                Wilcoxon.append("equal")
            Wilcoxon_p.append(p)

    return lists,normalty,normalty_p,Wilcoxon,Wilcoxon_p,taxa_list
"""
def two_groups(df):
    group1=[x*2 for x in range(int(df.shape[0]/2))]  #list1
    group2=[x*2-1 for x in range(1,int(df.shape[0]/2)+1)]  #list2
    return group1,group2
"""

# part 1: nonparametric test and output
def part1(df,group1,group1_name,group2,group2_name,paired,only_two_sided):
    # row1 and row2 belongs to group
    # group_1,group_2=two_groups(df)
    lists,normalty,normalty_p,Wilcoxon,Wilcoxon_p,taxa_list=\
    normalty_Wilcoxon_dual_group_test(df,group1,group2,paired,only_two_sided)

    output=[taxa_list,normalty,normalty_p,Wilcoxon,Wilcoxon_p]
    output_df=pd.DataFrame(output)
    output_df.index=["taxo","normalty","normalty_p","Wilcoxon","Wilcoxon_p"]
    if paired==True and only_two_sided!=True:
        name=str(filename+"_"+group1_name+"_vs_"+group2_name+".u-test.paired.csv")
    elif paired!=True and only_two_sided==True:
        name=str(filename+"_"+group1_name+"_vs_"+group2_name+".u-test.only_two_sided.csv")
    elif paired==True and only_two_sided==True:
        name=str(filename+"_"+group1_name+"_vs_"+group2_name+".u-test.paired.only_two_sided.csv")
    else:
        name=str(filename+"_"+group1_name+"_vs_"+group2_name+".u-test.csv")

    output_df.to_csv(name,header=None,sep=",")


# part 2: give the genes list which significant changed
# extract "increase", "decrease", or "equal" genes lists
def part2(df,target,group1,group1_name,group2,group2_name,groups_simple,groups,paired,only_two_sided):
    
    reordered_groups=[groups[x] for x in group1+group2]
    
    lists,normalty,normalty_p,Wilcoxon,Wilcoxon_p,taxa_list=\
    normalty_Wilcoxon_dual_group_test(df,group1,group2,paired,only_two_sided)
    abun_list=lists.copy() #[[[0.1,0.4,0.5...],[0.1,0.4,0.5...]],...]
    Wilcoxon_list=Wilcoxon.copy() #["increase","decrease"...]
    l2s=True    # large to small
    ceil_change=1

    targ_taxo_list=[]
    targ_abun_list=[]

    if len(abun_list)!=len(Wilcoxon_list):
        print("Error! length mismatch.")
        sys.exit()

    taxa_num=len(abun_list)

    for i in range(taxa_num):
        if Wilcoxon_list[i]==target:
            targ_taxo_list.append(taxa_list[i])
            targ_abun_list.append(abun_list[i])

    #if only_two_sided==False:
    if len(targ_abun_list)==0:
        print("no ",target," taxa. Exit")
        sys.exit()


    if l2s==True:
        # average relative change of certain taxa between after and before
        # (group2_abun_sum-group1_abun_sum)/group1_abun_sum
        # if group1_abun_sum==0, set the denominator as the minimum in group2_abun
        change_list=[]
        for x in targ_abun_list:
            try:
                relative_change=((sum(x[1])-sum(x[0]))/sum(x[0]))
                change_list.append(relative_change)
            except:
                # in case all abundance in denominator are zero
                denominator=min([value for value in x[1] if value!=0 ])
                relative_change=(sum(x[1])-sum(x[0]))/denominator
                change_list.append(relative_change)

        # targ_taxo_list: ['Bacteroides ilei']
        # targ_abun_list: the second item in taxo_abun_change
        # taxo_abun_change:
        # [('Bacteroides ilei', [[0.0, 0.003861003861003865, 0.002976190476190477, 0.0, 0.0, 0.0027752081406105405, 0.005994005994006001, 0.002450980392156865], [0.0, 0.010261194029850769, 0.00911577028258888, 0.00944206008583692, 0.012158054711246176, 0.010563380281690168, 0.0016260162601626016, 0.005016722408026762]], 2.2221268810078496)]
        taxo_abun_change=list(zip(targ_taxo_list,targ_abun_list,change_list))
        taxo_abun_change_l2s=taxo_abun_change.copy()
        # sort accordingly to the change
        taxo_abun_change_l2s.sort(reverse=True,key=lambda x:x[2])
        # max_change=max([x[2] for x in taxo_abun_change_l2s])
        # 10 % more than max average change for infinite change
        # ceil_change=max_change*1.1  # v1

    # output
    # Y-axis is sample ID. X-axis is taxa
    # file 1. the relative abundance
    # file 2. the relative change (not average change)
    # file 3. the average change

    # file 1
    ## relative abundance
    ## extract the abundance of each sample from taxo_abun_change_l2s
    ## in fact, the abundace in taxo_abun_change_l2s is the same as the original input table
    file1_y=reordered_groups.copy()
    file1_x=[x[0] for x in taxo_abun_change_l2s]
    ## x[1][0] are rel_abun of group1, x[1][1] are rel_abun of group 2

    file1_bulk=pd.DataFrame([x[1][0]+x[1][1] for x in taxo_abun_change_l2s]).T
    file1_df=file1_bulk.copy()
    file1_df.columns=file1_x
    file1_df.index=file1_y
    try:
        file1_df=file1_df.reindex(groups)
    except:
        return file1_df
    
    # file 2
    ## relative change
    # file2_y: ['305A', '396', '473B', '705E', '7108E', '711E', '778E', '782']
    ## if the increase change is more than one, set it as ceil_change
    ## the decrease change will never lower than -1
    if groups_simple !=None:
        file2_y=groups_simple.copy()
        file2_x=file1_x.copy()
        file2_bulk_list=[]
        for i in taxo_abun_change_l2s:
            relative_change_one_taxo=[]
            for j in range(len(file2_y)):
                try: # abs for negative
                    ## if the increase change is more than one, set it as ceil_change
                    ## the decrease change will never lower than -1
                    relative_change=min([ceil_change,(i[1][1][j]-i[1][0][j])/abs(i[1][0][j])]) #v3
                except:
                    print("float division by zero in ",i[0]," :")
                    if i[1][1][j]==0:
                        print("both numerator and denominator are zero.")
                        relative_change=0
                    elif  i[1][1][j]!=0:
                        print("use ceil change (10% more than max change)")
                        relative_change=ceil_change
                relative_change_one_taxo.append(relative_change)
            file2_bulk_list.append(relative_change_one_taxo)
    
        file2_bulk=pd.DataFrame(file2_bulk_list).T
        file2_df=file2_bulk.copy()
        file2_df.columns=file2_x
        file2_df.index=file2_y


    # file 3
    ## average change
    # (group2_abun_sum-group1_abun_sum)/group1_abun_sum
    # if group1_abun_sum==0, set the denominator as the minimum in group2_abun
    # not using ceil_change
    file3_y=["average_relative_change"]
    file3_x=file1_x.copy()
    """
    if target=="increase":
        file3_bulk=pd.DataFrame([x[2] for x in taxo_abun_change_l2s]).T
    elif target=="decrease":
        file3_bulk=pd.DataFrame([-x[2] for x in taxo_abun_change_l2s]).T
    """
    file3_bulk=pd.DataFrame([x[2] for x in taxo_abun_change_l2s]).T
    file3_df=file3_bulk.copy()
    file3_df.columns=file3_x
    file3_df.index=file3_y

    # write
    if only_two_sided==False:
        file1_name=str(filename+"_relative_abun_"+target+"_"+group1_name+"_vs_"+group2_name+".csv")
        file2_name=str(filename+"_relative_change_"+target+"_"+group1_name+"_vs_"+group2_name+".csv")
        file3_name=str(filename+"_average_change_"+target+".csv")

        file1_df.to_csv(file1_name,sep=',',index=True)
        if groups_simple !=None:
            file2_df.to_csv(file2_name,sep=',',index=True)
        file3_df.to_csv(file3_name,sep=',')
    elif only_two_sided==True:
        file1_name=str(filename+"_relative_abun_"+target+"_"+group1_name+"_vs_"+group2_name+".csv")
        file2_name=str(filename+"_relative_change_"+target+"_"+group1_name+"_vs_"+group2_name+".csv")
        file3_name=str(filename+"_average_change_"+target+".csv")

        file1_df.to_csv(file1_name,sep=',',index=True)
        if groups_simple !=None:
            file2_df.to_csv(file2_name,sep=',',index=True)
        file3_df.to_csv(file3_name,sep=',')



if  __name__=='__main__':
    import os
    if not os.path.exists('../taxa_abun/utest/'):
        os.makedirs('../taxa_abun/utest/')
        
    for i in range(1,9):
        assembly=sys.argv[1]
        group1=sys.argv[2]
        group1_name=sys.argv[3]
        group2=sys.argv[4]
        group2_name=sys.argv[5]
        
        dp=str(r"taxa_abun\rel_abun\sample18_abun_top.0."+str(i)+".rmU.euk.hvd.csv")
        
        basename=os.path.basename(dp)
        filename=os.path.join("../taxa_abun/utest/",basename)
        
        df=pd.read_csv(dp,sep=",",index_col=0)
        
        groups=list(df.index)
        #before
        ## paired=True
        group1=[6,7,8,9,10,11]
        group1_name="horse"
        #after
        group2=[0,1,2,3,4,5]
        group2_name="donkey"
        
        df.index=range(df.shape[0])
        only_two_sided=True
        paired=False
        ## groups_simple is None, when the samples are not paired
        #groups_simple=["305A","396","473B","705E","7108E","711E","778E","782"] 
        groups_simple=None

        ## group2 comparing to group1


        """
        #SC
        group1=[0,1,2,3,8,9,10,11]
        group1_name="SC"
        #TC
        group2=[4,5,6,7,12,13,14,15]
        group2_name="TC"
        """

        #cow
        """
        group1=[0,1,2,3,4,5]
        group1_name="cow"
        #heifer
        group2=[6,7,8,9,10,11,12,13,14,15]
        group2_name="heifer"
        """


        ## if paired==True, use scipy.stats.wilcoxon paired test
        ## else use stats.mannwhitneyu

        group1_group2_names=[group1_name]*len(group1)+[group2_name]*len(group2)

        if only_two_sided==True:
            for target in ["equal","unequal"]:
                part1(df,group1,group1_name,group2,group2_name,paired,only_two_sided)
                part2(df,target,group1,group1_name,group2,group2_name,groups_simple,groups,paired,only_two_sided)
        if only_two_sided!=True:
            for target in ["increase","decrease"]:
                part1(df,group1,group1_name,group2,group2_name,paired,only_two_sided)
                part2(df,target,group1,group1_name,group2,group2_name,groups_simple,groups,paired,only_two_sided)

