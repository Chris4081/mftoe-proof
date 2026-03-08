import pandas as pd

df = pd.read_csv(
    "Pantheon+SH0ES.dat",
    sep=r"\s+",
    comment="#"
)

print(df.columns)

df2 = df[["zHD", "MU_SH0ES", "MU_SH0ES_ERR_DIAG"]]

df2.columns = ["z", "mu", "mu_err"]

df2.to_csv("pantheon_plus.csv", index=False)

print("Saved pantheon_plus.csv")