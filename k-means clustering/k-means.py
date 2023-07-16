import tkinter as tk
from tkinter import filedialog
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import math

root= tk.Tk()


canvas1 = tk.Canvas(root, width = 400, height = 300,  relief = 'raised')
canvas1.pack()

label1 = tk.Label(root, text='k-Means Clustering')
label1.config(font=('helvetica', 14))
canvas1.create_window(200, 25, window=label1)

label2 = tk.Label(root, text='Type Number of Clusters:')
label2.config(font=('helvetica', 8))
canvas1.create_window(200, 120, window=label2)

entry1 = tk.Entry (root) 
canvas1.create_window(200, 140, window=entry1)

def getExcel ():
    
    global df
    import_file_path = filedialog.askopenfilename()
    read_file = pd.read_excel (import_file_path)
    df = DataFrame(read_file,columns=['x','y'])  
    
browseButtonExcel = tk.Button(text=" Import Excel File ", command=getExcel, bg='green', fg='white', font=('helvetica', 10, 'bold'))
canvas1.create_window(200, 70, window=browseButtonExcel)

def getKMeans ():
    global df
    global numberOfClusters
    numberOfClusters = int(entry1.get())

    kmeans = KMeans(n_clusters=numberOfClusters, init="k-means++").fit(df)
    centroids = kmeans.cluster_centers_
    
    ar = kmeans.labels_
    print (kmeans.inertia_)
    
    i = 0
    file2 = str(len(ar)) + "Sensors - " + str(numberOfClusters) + "cluster - size of clusters.txt"
    f2 = open(file2, "w")
    while i < numberOfClusters :
        file = str(len(ar)) + "Sensors - " + str(numberOfClusters) + "cluster - " + str(i + 1) + ".tsp"
        f = open(file, "w")
        j = 0 
        k = 1
        s = "NAME : " + file + "\n" + "COMMENT : 51-city problem (Christofides/Eilon)\nTYPE : TSP\nDIMENSION : 51\nEDGE_WEIGHT_TYPE : EUC_2D\nNODE_COORD_SECTION\n"
        f.write(s)
        while j < len(ar) :
            if(ar[j] == i) :
                list = df.iloc[j].dropna().tolist()
                a = list[0] * 100
                b = list[1] * 100
                p = str(k) + " " + str(round(a)) + " " + str(round(b)) + "\n"
                f.write(p)
                k += 1
            j += 1
        p2 = "Cluster " + str(i + 1) + " : " + str(k - 1) + "\n"
        f2.write(p2)
        i += 1
        f.write("EOF")
      
    label3 = tk.Label(root, text= centroids)
    file3 = str(len(ar)) + "Sensors - " + str(numberOfClusters) + "cluster - Centroids of clusters.txt"
    f3 = open(file3, "w")
    string = str(centroids)
    f3.write(string)
    f3.write("\nEOF")
    canvas1.create_window(500, 500, window=label3)
    
    figure1 = plt.Figure(figsize=(5,5), dpi=100)
    ax1 = figure1.add_subplot()
    scatter = ax1.scatter(df['x'], df['y'], c= kmeans.labels_.astype(float), s=50, alpha=0.5)
    ax1.scatter(centroids[:, 0], centroids[:, 1], c='red', s=50)
    scatter1 = FigureCanvasTkAgg(figure1, root) 
    scatter1.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH)
    
    
processButton = tk.Button(text=' Process k-Means ', command=getKMeans, bg='brown', fg='white', font=('helvetica', 10, 'bold'))
canvas1.create_window(200, 170, window=processButton)



root.mainloop()
