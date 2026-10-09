

import string

import pandas as pd
import pathlib as path
import sys
repo=path.Path(__file__).resolve().parent.parent
results_path=repo/'results'/'polar_50c.csv'
results_check=results_path.exists()
n_last=500
rel_tol=0.001
abs_floor=1e-6  



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


def convergence_check(df):
    if len(df) < n_last:
        print(f"Dataframe has less than {n_last} rows. Convergence check cannot be performed.")
        return False
    all_converged=True
    # Check for convergence in the last n_last rows of the dataframe with Cl, Cd, and Cm columns
    last_rows=df.tail(n_last)
    for col in ['Cl', 'Cd', 'Cm']:
        mean=last_rows[col].mean()
        maximum=last_rows[col].max()
        minimum=last_rows[col].min()
        spread=maximum-minimum
        print(f"{col}: Mean: {mean}, Spread: {spread}")
        allowance=max(rel_tol*abs(mean), abs_floor)
        spread_percentage=(spread/allowance)*100
        print(f"Spread as a percentage of allowance: {spread_percentage:.2f}%")
        if spread_percentage >=100:
            print(f"{col} has not converged. Spread {spread} exceeds allowance {allowance}.")
            all_converged=False
    return all_converged
    


summary_table=[]

AoA=[0,2,4,6,8,10,12]
for angle in AoA:
    data_path=repo/'case'/'results'/f'AoA{angle}'/'postProcessing'/'forceCoeffs1'/'0'/'forceCoeffs.dat'
    if not data_path.exists():
        sys.exit(f"Data file not found at {data_path}")
    data_output=reader_function(data_path)
    converged=convergence_check(data_output)

    print(f"Convergence check for AoA {angle}: {'PASS' if converged else 'FAIL'}")
    result=file_read[file_read['alpha_deg']==angle]
    final_result=result.iloc[0]
    compare_Cl=abs(final_result['CL']-data_output.iloc[-1]['Cl'])
    compare_Cd=abs(final_result['CD']-data_output.iloc[-1]['Cd'])
    compare_Cm=abs(final_result['Cm_c4']-data_output.iloc[-1]['Cm'])
    print(f"Comparison of coefficients for AoA {angle}:")
    print(f"  Cl: {compare_Cl:.1e}")
    print(f"  Cd: {compare_Cd:.1e}")
    print(f"  Cm: {compare_Cm:.1e}")
    max_diff=max(compare_Cl, compare_Cd, compare_Cm)
    matches_polarCsv=max_diff<=5e-7
    print(f"  Max difference: {max_diff:.1e}")
    summary_table.append({
        'alpha': angle,
        'cl': data_output.iloc[-1]['Cl'],
        'cd': data_output.iloc[-1]['Cd'],
        'cm': data_output.iloc[-1]['Cm'],
        'max_diff': max_diff,
        'matches_polarCsv': matches_polarCsv,
        'converged': converged
    })
    if matches_polarCsv:
        print("Comparison check: PASS")
    else:
        print("Comparison check: FAIL")

final_results=pd.DataFrame(summary_table).to_string(index=False)

print(final_results)








