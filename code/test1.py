import pandas as pd

F = pd.read_excel(r"D:\workspace\SHUWEICUP2544531_2025.11.14\ProblemA\timeAndForce.xlsx", usecols='B')

print(F.loc[999])