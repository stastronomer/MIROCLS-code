import split_spitzer_data
import cluster_searcher_2


print("\n \n----Mid-IR OPTICS CLuster Search (MIROCLS)---- \n \n Masterarbeit @Akash Gupta @European Southern Observatory @Ludwig Maximilians University, Munich\n \n")
print("1. Press and enter S to run the splitting algorithm")
print("2. Press and enter C to run the Cluster Search algorithm using sklearn.optics")
print("3. Press and enter V to run the visualizer code to plot the CMDs and scatter plot")
print("4. Press and enter A to run the whole package step by step (achtung)")

k = input()

if k=='S':
    split_spitzer_data.split()
elif k=='C':
    split_spitzer_data.split()
elif k=='V':
    import visualizer.py

elif k == 'A':
    print("Starting splitting...\n")
    split_spitzer_data.split()
    print("\n \n Spltting done.")

    print("\n Clustering on split data \n")
    split_spitzer_data.split()
    print("\n Clustering done... \n")

    print("\n Starting visualizer code... \n")
    import visualizer.py
    print("\n Plots done. \n \n")

else:
    print("Please select appropriate prompt.")

print("Thanks for using MIROCLS code.")