from __main__ import *
#import pandas as pd




file=data_path+"alphas_table.csv"
df_selection_params=pd.read_csv(file,sep=",")#,index_col=0)
display(df_selection_params)

file2=data_path+"parameters_table.csv"
df_growth_params=pd.read_csv(file2)
display(df_growth_params)




NpBGT=19
#NpBGT=35
print("Number of plasmids",NpBGT)

amps=list(df_selection_params["Amp"])
kappaN=df_selection_params["κ_n"][0]
alphas=list(df_selection_params["⍺"])
print("kappaN",kappaN)
print("alphas",alphas)

r=float(df_growth_params[df_growth_params["Parameter"]=="r"]["Measured value"][0])
print("r: ",r)
ρ=float(list(df_growth_params[df_growth_params["Parameter"]=="ρ"]["Measured value"])[0])

ρN=ρ/NpBGT
print("ρ:",ρ)
print("ρN:",ρN)
σ=float(list(df_growth_params[df_growth_params["Parameter"]=="σ"]["Measured value"])[0])
print("σ",σ)

pBGT_cost=kappaN/NpBGT
print("Cost per plasmid:",pBGT_cost)

mu_n=float(df_growth_params[df_growth_params["Parameter"]=="μ_n"]["Estimated value"])



