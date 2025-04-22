import pandas as pd
import matplotlib.pyplot as plt

# Load CSV and clean column names
file_name = input("Enter the name of the CSV file (with .csv extension): ")
proteins = input("Enter the proteins belong to chromosomal/plasmid :")
df = pd.read_csv(file_name)
df.columns = [col.strip() for col in df.columns]

# Count mutations
mutation_counts = df["Number of Mutations"]
mutated = (mutation_counts > 0).sum()
non_mutated = (mutation_counts == 0).sum()
total_proteins = len(df)

# Pie chart data
labels = ["Mutated Proteins", "Non-Mutated Proteins"]
sizes = [mutated, non_mutated]
colors = ["#FF6F61", "#6BAED6"]
explode = (0.1, 0)  # Detach mutated slice

# Plot pie chart
plt.figure(figsize=(8, 8))
plt.pie(sizes, labels=labels, colors=colors, explode=explode,
        autopct='%1.1f%%', startangle=140, shadow=True,
        textprops={'fontsize': 14})
plt.title(f"Mutation Status Among {total_proteins} {proteins} Proteins\n"
          f"Mutated: {mutated} | Non-Mutated: {non_mutated}",
          fontsize=15, fontweight='bold')
plt.tight_layout()
plt.show()