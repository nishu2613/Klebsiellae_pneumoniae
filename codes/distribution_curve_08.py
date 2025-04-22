import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Step 1: Load the CSV file
file_path = "distribution_curve_plasmid.csv"  # Change this if your file name is different
df = pd.read_csv(file_path)

# Step 2: Normalize column names
df.columns = [col.strip().lower() for col in df.columns]

# Step 3: Rename columns for consistency (optional, based on your current column names)
df = df.rename(columns={
    'id': 'protein_id',
    'length': 'length',
    'mutation_status': 'mutation_status'
})

# Step 4: Clean 'length' column and drop invalid entries
df['length'] = pd.to_numeric(df['length'], errors='coerce')
df = df.dropna(subset=['length'])

# Step 5: Normalize mutation_status values (make lowercase, remove spaces)
df['mutation_status'] = df['mutation_status'].str.strip().str.lower()

# Step 6: Plot KDE (Density) curve for mutation vs no mutation
plt.figure(figsize=(10, 6))

# KDE for mutated proteins
sns.kdeplot(
    data=df[df["mutation_status"] == "mutation"],
    x="length", label="Mutated", fill=True,
    color="red", alpha=0.6
)

# KDE for non-mutated proteins
sns.kdeplot(
    data=df[df["mutation_status"] == "no mutation"],
    x="length", label="Non-mutated", fill=True,
    color="green", alpha=0.6
)

# Step 7: Plot styling
plt.title("Protein Length Distribution: Mutated vs Non-mutated")
plt.xlabel("Protein Length (amino acids)")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.show()