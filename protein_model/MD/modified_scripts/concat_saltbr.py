import os
import re
import pandas as pd

df_list = []
for fname in os.listdir("."):
    m1 = re.findall(r"saltbr-([A-Z]{3}\d+)_chain([AB])-([A-Z]{3}\d+)_chain([AB])\.dat", fname)
    m2 = re.findall(r"saltbr-([A-Z]{3}\d+)-([A-Z]{3}\d+)\.dat", fname)
    if m1:
        a1, c1, a2, c2 = m1[0]
        df = pd.read_csv(fname, header=None, sep=" ")
        df.columns = ["frame", "distance"]
        df["recipient"] = a1
        df["recipient_chain"] = c1
        df["donor"] = a2
        df["donor_chain"] = c2
        df_list.append(df)
    elif m2:
        a1, a2 = m2[0]
        c1 = "A" if int(a1[3:]) < 598 else "B"
        c2 = "A" if int(a2[3:]) < 598 else "B"
        df = pd.read_csv(fname, header=None, sep=" ")
        df.columns = ["frame", "distance"]
        df["recipient"] = a1
        df["recipient_chain"] = c1
        df["donor"] = a2
        df["donor_chain"] = c2
        df_list.append(df)

pd.concat(df_list).to_csv("saltbr-all.csv", index=False)