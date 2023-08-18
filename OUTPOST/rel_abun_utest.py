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
,Clostridiales,bacterium,Verrucomicrobia,bacterium
donkey1,0.03258823,0.017380389
donkey2,0.035035599,0.016954362
donkey3,0.042927956,0.015255787

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
import os

def extract_lists(df,column,index1,index2):
    list1=df.loc[index1,column].tolist()
    list2=df.loc[index2,column].tolist()
    return [list1,list2]

def normalty_Wilcoxon_dual_group_test(df,group1,group2,paired=False,only_two_sided=True):
    # lists of abundance, list1 and list2 responds to group1 and group2
    abun_list=[] #[[list1,list2],[]....]
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
        [list1,list2]=extract_lists(df,col,group1,group2)
        list12=list1+list2
        abun_list.append([list1,list2])

        alpha=1e-3
        k, p = stats.kstest(list12,'norm')
        if p<alpha:
            normalty.append("N")
            normalty_p.append(p)
        else:
            normalty.append("Y")
            normalty_p.append(p)

        alpha=0.05
        if only_two_sided == False:
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

    return abun_list,normalty,normalty_p,Wilcoxon,Wilcoxon_p,taxa_list

# part 1: nonparametric test and output
def part1(df,group1,group1_name,group2,group2_name,paired,only_two_sided):
    # row1 and row2 belongs to group
    # group_1,group_2=two_groups(df)
    _,normalty,normalty_p,Wilcoxon,Wilcoxon_p,taxa_list = \
    normalty_Wilcoxon_dual_group_test(df,group1,group2,paired,only_two_sided)

    output=[taxa_list,normalty,normalty_p,Wilcoxon,Wilcoxon_p]
    output_df=pd.DataFrame(output)
    output_df.index=["taxo","normalty","normalty_p","Wilcoxon","Wilcoxon_p"]
    if paired==True and only_two_sided!=True:
        name=str(prefix+".u-test.paired.csv")
    elif paired!=True and only_two_sided==True:
        name=str(prefix+".u-test.two_sided.csv")
    elif paired==True and only_two_sided==True:
        name=str(prefix+".u-test.paired.two_sided.csv")
    else:
        name=str(prefix+".u-test.csv")

    output_df.to_csv(name,header=None,sep=",")


# part 2: give the genes list which significant changed
# extract ("increase","decrease"), or ("equal","unequal") genes lists
def part2(df,target,group1,group1_name,group2,group2_name,subgroup,groups,paired,only_two_sided):

    # reordered_groups = [groups[x] for x in group1+group2]

    abun_list,normalty,normalty_p,Wilcoxon,Wilcoxon_p,taxa_list =\
    normalty_Wilcoxon_dual_group_test(df,group1,group2,paired,only_two_sided)
    Wilcoxon_list = Wilcoxon.copy() #["equal","unequal"...]
    ceil_change = 1

    targ_taxo_list = []
    targ_abun_list = []

    assert len(abun_list) == len(Wilcoxon_list),"OUTPOST: rel_abun_utest: length mismatch. exit."

    # print(f"Wilcoxon_list:{Wilcoxon_list}")
    # print(f"group1_name:{group1_name}")
    for i in range(len(abun_list)):
        if Wilcoxon_list[i] == target:
            targ_taxo_list.append(taxa_list[i])
            targ_abun_list.append(abun_list[i])
    # print(f"targ_taxo_list:{targ_taxo_list}")
    # print(f"targ_abun_list:{targ_abun_list}")

    # write
    file1_name = str(prefix+".rel_abun."+target+".csv")
    file2_name = str(prefix+".rel_change."+target+".csv")
    file3_name = str(prefix+".ave_change."+target+".csv")

    if len(targ_abun_list) == 0:
        print("OUTPOST: rel_abun_utest: no ",target," taxa.")
        os.system("touch {}".format(file1_name))
        os.system("touch {}".format(file2_name))
        os.system("touch {}".format(file3_name))
        return

    # order to large to small
    # average relative change of certain taxa comparison
    # (group2_abun_sum-group1_abun_sum)/group1_abun_sum
    # if group1_abun_sum==0, set the denominator as the minimum in group2_abun
    change_list = []
    for x in targ_abun_list:
        try:
            rel_change = ((sum(x[1])-sum(x[0]))/sum(x[0]))
            change_list.append(rel_change)
        except:
            # in case all abundance in denominator are zero
            # print(f"X: {x}")
            # print(f"target: {target}")
            try:
                denominator = min([value for value in x[1] if value!=0 ])
                rel_change = (sum(x[1])-sum(x[0]))/denominator
                change_list.append(rel_change)
            except:
                change_list.append(0)

    # targ_taxo_list: ['Bacteroides ilei']
    # targ_abun_list: the second item in taxo_abun_change
    # taxo_abun_change:
    # [('Bacteroides ilei', [[0.0, 0.003861003861003865, 0.002976190476190477, 0.0, 0.0, 0.0027752081406105405, 0.005994005994006001, 0.002450980392156865], [0.0, 0.010261194029850769, 0.00911577028258888, 0.00944206008583692, 0.012158054711246176, 0.010563380281690168, 0.0016260162601626016, 0.005016722408026762]], 2.2221268810078496)]
    assert len(targ_taxo_list) == len(targ_abun_list) == len(change_list),\
    "OUTPOST: rel_abun_utest: length mismatch. exit."

    taxo_abun_change = list(zip(targ_taxo_list,targ_abun_list,change_list))
    taxo_abun_change_l2s = taxo_abun_change.copy()
    # sort accordingly to the change
    taxo_abun_change_l2s.sort(reverse = True,key = lambda x:x[2])
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
    # file1_y = reordered_groups.copy()
    file1_x = [x[0] for x in taxo_abun_change_l2s]

    ## x[1][0] are rel_abun of group1, x[1][1] are rel_abun of group 2
    file1_df = pd.DataFrame([x[1][0]+x[1][1] for x in taxo_abun_change_l2s]).T
    file1_df.columns = file1_x
    file1_df.index = [groups[x] for x in group1+group2]
    # file1_df = file1_df.reindex(groups)
    # file 2
    ## relative change
    # file2_y: ['305A', '396', '473B', '705E', '7108E', '711E', '778E', '782']
    ## if the increase change is more than one, set it as ceil_change
    ## the decrease change will never lower than -1
    if subgroup !=None:
        file2_y = subgroup.copy()
        file2_x = file1_x
        file2_bulk_list=[]
        for i in taxo_abun_change_l2s:
            rel_change_one_taxo=[]
            for j in range(len(file2_y)):
                try: # abs for negative
                    ## if the increase change is more than one, set it as ceil_change
                    ## the decrease change will never lower than -1
                    rel_change=min([ceil_change,(i[1][1][j]-i[1][0][j])/abs(i[1][0][j])]) #v3
                except:
                    print("float division by zero in ",i[0]," :")
                    if i[1][1][j]==0:
                        print("both numerator and denominator are zero.")
                        rel_change=0
                    elif  i[1][1][j]!=0:
                        print("use ceil change (10% more than max change)")
                        rel_change=ceil_change
                rel_change_one_taxo.append(rel_change)
            file2_bulk_list.append(rel_change_one_taxo)

        file2_bulk=pd.DataFrame(file2_bulk_list).T
        file2_df=file2_bulk.copy()
        file2_df.columns=file2_x
        file2_df.index=file2_y


    # file 3
    ## average change
    # (group2_abun_sum-group1_abun_sum)/group1_abun_sum
    # if group1_abun_sum==0, set the denominator as the minimum in group2_abun
    # not using ceil_change
    file3_y=["average_rel_change"]
    file3_x = file1_x
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



    file1_df.to_csv(file1_name,sep=',',index=True)
    if subgroup != None:
        file2_df.to_csv(file2_name,sep=',',index=True)
    file3_df.to_csv(file3_name,sep=',')


if  __name__ == '__main__':
    try:
        dp_list = sys.argv[1].split(',')
        group1 = [int(x) for x in sys.argv[2].split(',')]
        group1_name = sys.argv[3]
        group2 = [int(x) for x in sys.argv[4].split(',')]
        group2_name = sys.argv[5]
        prefix_list = sys.argv[6].split(',')
        
        try:
            paired = eval(sys.argv[7])
            only_two_sided = eval(sys.argv[8])
        except:
            sys.exit("OUTPOST: rel_abun_utest: paired/two-sided parameter wrong. exit.")
        assert paired in [True, False] and only_two_sided in [True,False], "OUTPOST: rel_abun_utest: paired/two-sided parameter wrong. exit."

        for dp,prefix in zip(dp_list,prefix_list):
            df = pd.read_csv(dp,sep = ",",index_col = 0)

            groups = list(df.index)

            df.index = range(df.shape[0])
            #only_two_sided = True
            ## subgroup is None, when the samples are not paired
            #subgroup=["305A","396","473B","705E","7108E","711E","778E","782"]
            subgroup = None

            ## group2 comparing to group1
            ## if paired==True, use scipy.stats.wilcoxon paired test
            ## else use stats.mannwhitneyu

            if only_two_sided == True:
                for target in ["equal","unequal"]:
                    part1(df,group1,group1_name,group2,group2_name,paired,only_two_sided)
                    part2(df,target,group1,group1_name,group2,group2_name,subgroup,groups,paired,only_two_sided)
            else:
                for target in ["increase","decrease"]:
                    part1(df,group1,group1_name,group2,group2_name,paired,only_two_sided)
                    part2(df,target,group1,group1_name,group2,group2_name,subgroup,groups,paired,only_two_sided)
    except Exception as e:
        import traceback
        error_log = sys.argv[9]
        os.system("touch " + error_log)
        print(f"OUTPOST: {e}")
        traceback.print_exc()
