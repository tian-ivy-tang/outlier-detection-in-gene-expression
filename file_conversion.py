import pandas as pd

file = pd.read_csv("KremerNBaderSmall.tsv", sep="\t")
print(file.head())
print(file.shape)

def add_one(x):
    try:
        return x + 1
    except TypeError:
        return x

# Apply the function to each cell using map
file = file.map(add_one)
print(file.head())
print(file.shape)
file.to_csv('KremerNBaderSmall.csv', sep='\t')