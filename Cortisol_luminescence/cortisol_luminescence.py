import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from scipy import stats
import pandas as pd
import csv

# load file
df = pd.read_excel(r'D:\Olfaction\cortisol_luminescence\internal_cortisol_Sleep_deprivation\cortisol_measured.xlsx')
df = df.dropna()
df.Treatment.unique()
print(df.Treatment.unique())
df = df[(df.Treatment=='control ')| (df.Treatment=='light pulse ')|(df.Treatment=='continous')|
       (df.Treatment=='control2')|(df.Treatment=='voltage5V_15min')]

print(f"control mean: {df[(df.Treatment=='control ')].mean()} control std: {df[(df.Treatment=='control ')].std()}")
print(f"continous mean: {df[(df.Treatment=='continous')].mean()} continous std: {df[(df.Treatment=='continous')].std()}")
print(f"light pulse mean: {df[(df.Treatment=='light pulse ')].mean()} light pulse std: {df[(df.Treatment=='light pulse ')].std()}")
print(f"control 2 mean: {df[(df.Treatment=='control2')].mean()} control 2 std: {df[(df.Treatment=='control2')].std()}")
print(f"voltage mean: {df[(df.Treatment=='voltage5V_15min')].mean()} voltage std: {df[(df.Treatment=='voltage5V_15min')].std()}")

# load second file
# focus on free swimming 30 minute treatment
df_cort = pd.read_excel(r'D:\Olfaction\cortisol_luminescence\internal_cortisol_cortisol_treatment\interpolated_readout_may2023.xlsx')
df_cort = df_cort.dropna()
df_cort.Treatment.unique()
print(df_cort.Treatment.unique())

# create df just for subgroups
df_int = df_cort[['Treatment','cortisol (ug/dl)','cortisol (ng/ml)']].copy()
df_fs30 = df_int[df_int['Treatment'].str.startswith('FS30mi')]
print(df_fs30)

#concatenate data frames
df3 = df.append(df_fs30, ignore_index=True)
print(df3)

## statistics
#define samples
group1 = df[df['Treatment']=='control ']
group2 = df[df['Treatment']=='light pulse ']
#perform independent two sample t-test
p1=stats.ttest_ind(group1['cortisol (ng/ml)'], group2['cortisol (ng/ml)'],alternative ='less')
print("lightpulse",p1)

group1 = df[df['Treatment']=='control ']
group2 = df[df['Treatment']=='continous']
#perform independent two sample t-test
p2=stats.ttest_ind(group1['cortisol (ng/ml)'], group2['cortisol (ng/ml)'], alternative ='less')
print("continuous",p2)

group1 = df[df['Treatment']=='control2']
group2 = df[df['Treatment']=='voltage5V_15min']
#perform independent two sample t-test
p3=stats.ttest_ind(group1['cortisol (ng/ml)'], group2['cortisol (ng/ml)'],alternative ='less')
print("voltage",p3)