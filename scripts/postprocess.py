from modulefinder import test

import pandas as pd
import pathlib as path
import sys
repo=path.Path(__file__).resolve().parent.parent
results_path=repo/'results'/'polar_50c.csv'


results_check=results_path.exists()




if not results_check:
    sys.exit(f"Results file not found at {results_path}")

#----------------------------------------------

def reader_function(file_path):
    data_lines=[]
    with open(file_path) as f:
        for l in f: 
            if l.startswith("#"):
                data_lines.append(l)
    header=data_lines[-1]
    stripped_header=header.lstrip("#")
    cols=stripped_header.split()
    df=pd.read_csv(file_path,comment="#", sep=r"\s+", names=cols)
    return df

file_read=pd.read_csv(results_path)

AoA=[0,2,4,6,8,10,12]
for angle in AoA:
    data_path=repo/'case'/'results'/f'AoA{angle}'/'postProcessing'/'forceCoeffs1'/'0'/'forceCoeffs.dat'
    if not data_path.exists():
        sys.exit(f"Data file not found at {data_path}")
    data_output=reader_function(data_path)
    result=file_read[file_read['alpha_deg']==angle]
    final_result=result.iloc[0]
    compare_Cl=abs(final_result['CL']-data_output.iloc[-1]['Cl'])
    compare_Cd=abs(final_result['CD']-data_output.iloc[-1]['Cd'])
    compare_Cm=abs(final_result['Cm_c4']-data_output.iloc[-1]['Cm'])
    print(f"Comparison for AoA {angle}:")
    print(f"  Cl: {compare_Cl:.1e}")
    print(f"  Cd: {compare_Cd:.1e}")
    print(f"  Cm: {compare_Cm:.1e}")
    max_diff=max(compare_Cl, compare_Cd, compare_Cm)
    print(f"  Max difference: {max_diff:.1e}")
    if max_diff <=5e-7:
        print("Result: PASS")
    else:
        print("Result: FAIL")








#Reader function to read OpenFOAM data files

#data_lines=[]
#with open(data_path) as f:
#    for l in f: 
#        if l.startswith("#"):
#            data_lines.append(l)

#header=data_lines[-1]
#stripped_header=header.lstrip("#")
#cols=stripped_header.split()



#Output the dataframe and its properties
###print(df)
#print(df.shape)
#print(df.columns)
#print(df.tail(1))

#print(df['Time'].head(1))
#print(df['Time'].tail(1))
#print(df.head(3))

