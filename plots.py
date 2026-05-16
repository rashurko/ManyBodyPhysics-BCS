import matplotlib.pyplot as plt
import numpy as np

def plotTotEnergy(dataList, labels, colors, markers, linestyles, title, markersize=0, name="default.png"):
    fig, ax = plt.subplots()
    for data, label, color, marker, linestyle in zip(dataList, labels, colors, markers, linestyles):
        ax.plot(data[:,0], data[:,1], color=color, marker=marker, ms=markersize, linestyle=linestyle, label=label)
    # Calculate the difference between energies where they have the same g values
    diff = []
    gDiff = []
    for i in range(len(dataList[0])):
        for j in range(len(dataList[1])):
            if dataList[0][i, 0] == dataList[1][j, 0]:
                diffE = dataList[0][i, 1] - dataList[1][j, 1]
                gDiff.append(dataList[0][i, 0])
                diff.append(diffE)
    ax.plot(gDiff, diff, color="red", linestyle="dashdot", label=r"$\Delta E_0$")
    ax.set_xlabel("g")
    ax.set_ylabel(r"$E_0$")
    ax.set_title(title)
    plt.legend()
    plt.savefig(name)
    plt.close()

def plotEp1(dataList, labels, colors, markers, linestyles, title, markersize=0, name="default.png"):
    fig, ax = plt.subplots()
    for data, label, color, marker, linestyle in zip(dataList, labels, colors, markers, linestyles):
        ax.plot(data[:,0], data[:,2], color=color, marker=marker, ms=markersize, linestyle=linestyle, label=label)
    diff = []
    gDiff = []
    for i in range(len(dataList[0])):
        for j in range(len(dataList[1])):
            if dataList[0][i, 0] == dataList[1][j, 0]:
                diffE = dataList[0][i, 2] - dataList[1][j, 2]
                gDiff.append(dataList[0][i, 0])
                diff.append(diffE)
    ax.plot(gDiff, diff, color="red", linestyle="dashdot", label=r"$\Delta E_{p1}$")
    ax.set_xlabel("g")
    ax.set_ylabel(r"$E_{p1}$")
    ax.set_title(title)
    plt.legend()
    plt.savefig(name)
    plt.close()

def plotEpair(dataList, labels, colors, markers, linestyles, title, markersize=0, name="default.png"):
    fig, ax = plt.subplots()
    for data, label, color, marker, linestyle in zip(dataList, labels, colors, markers, linestyles):
        ax.plot(data[:,0], data[:,3], color=color, marker=marker, ms=markersize, linestyle=linestyle, label=label)
    diff = []
    gDiff = []
    for i in range(len(dataList[0])):
        for j in range(len(dataList[1])):
            if dataList[0][i, 0] == dataList[1][j, 0]:
                diffE = dataList[0][i, 1] - dataList[1][j, 1]
                gDiff.append(dataList[0][i, 0])
                diff.append(diffE)
    ax.plot(gDiff, diff, color="red", linestyle="dashdot", label=r"$\Delta E_{pair}$")
    ax.set_xlabel("g")
    ax.set_ylabel(r"$E_{pair}$")
    ax.set_title(title)
    plt.legend()
    plt.savefig(name)
    plt.close()

def plotDuration(dataList, labels, colors, markers, linestyles, title, markersize=0, name="default.png"):
    fig, ax = plt.subplots()
    for data, label, color, marker, linestyle in zip(dataList, labels, colors, markers, linestyles):
        ax.plot(data[:,0], data[:,4], color=color, marker=marker, ms=markersize, linestyle=linestyle, label=label)
    ax.set_xlabel("g")
    ax.set_ylabel(r"Duration ($\mu$s)")
    ax.set_title(title)
    plt.legend()
    plt.savefig(name)
    plt.close()

def plotOccupationHistogram(dataList, labels, colors, title, name="default.png"):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for data, label, color in zip(dataList, labels, colors):
        for i in range(0, len(data)):
            if i == 0:
                ax.bar3d(np.arange(0, len(data[i,1:])), data[i,0], 0, dx=1, dy=0.01, dz=data[i,1:], color=color, alpha=0.5, label=label)
            else:
                ax.bar3d(np.arange(0, len(data[i,1:])), data[i,0], 0, dx=1, dy=0.01, dz=data[i,1:], color=color, alpha=0.5)
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(r"$g$")
    ax.set_zlabel(r"$\rho _\alpha$")
    ax.set_title(title)
    plt.legend()
    plt.savefig(name)
    plt.close()

def plotDiffOccupationHistogram(dataList, title, name="default.png"):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    rhoDiff = []
    gDiff = []
    for i in range(len(dataList[0])):
        for j in range(len(dataList[1])):
            if dataList[0][i, 0] == dataList[1][j, 0]:
                rhoDiff.append(np.abs(dataList[0][i, :] - dataList[1][j, :]))
                gDiff.append(dataList[0][i, 0])
    
    for i in range(len(rhoDiff)):
        if i == 0:
            ax.bar3d(np.arange(0, len(rhoDiff[i][1:])), gDiff[i], 0, dx=1, dy=0.01, dz=rhoDiff[i][1:], color="red", alpha=0.5, label=r"$|\Delta \rho _\alpha|$")
        else:
            ax.bar3d(np.arange(0, len(rhoDiff[i][1:])), gDiff[i], 0, dx=1, dy=0.01, dz=rhoDiff[i][1:], color="red", alpha=0.5)
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(r"$g$")
    ax.set_zlabel(r"$|\Delta \rho _\alpha|$")
    ax.set_title(title)
    # Set pov
    ax.view_init(elev=15, azim=-15)
    plt.legend()
    plt.savefig(name)
    plt.close()

# BCS data
bcs_data_N4 = np.loadtxt("Data/BCS/N4/data.txt", skiprows=1)
bcs_data_N12 = np.loadtxt("Data/BCS/N12/data.txt", skiprows=1)
bcs_occupation_N4 = np.loadtxt("Data/BCS/N4/occupation.txt", skiprows=1)
bcs_occupation_N12 = np.loadtxt("Data/BCS/N12/occupation.txt", skiprows=1)
# Exact data
exact_data_N4 = np.loadtxt("Data/Exact/N4/data.txt", skiprows=1)
exact_data_N12 = np.loadtxt("Data/Exact/N12/data.txt", skiprows=1)
exact_occupation_N4 = np.loadtxt("Data/Exact/N4/occupation.txt", skiprows=1)
exact_occupation_N12 = np.loadtxt("Data/Exact/N12/occupation.txt", skiprows=1)

plotTotEnergy([bcs_data_N4, exact_data_N4], ["BCS", "Exact"], ["black", "black"], [".", "x"], ["-", "--"], r"$E_0$ BCS vs Exact (N=4)", markersize=0, name="Plots/E0_N4")
plotTotEnergy([bcs_data_N12, exact_data_N12], ["BCS", "Exact"], ["black", "black"], [".", "x"], ["-", "--"], r"$E_0$ BCS vs Exact (N=12)", markersize=0, name="Plots/E0_N12")

plotEp1([bcs_data_N4, exact_data_N4], ["BCS", "Exact"], ["black", "black"], [".", "x"], ["-", "--"], r"$E_{p1}$ BCS vs Exact (N=4)", markersize=0, name="Plots/Ep1_N4")
plotEp1([bcs_data_N12, exact_data_N12], ["BCS", "Exact"], ["black", "black"], [".", "x"], ["-", "--"], r"$E_{p1}$ BCS vs Exact (N=12)", markersize=0, name="Plots/Ep1_N12")

plotEpair([bcs_data_N4, exact_data_N4], ["BCS", "Exact"], ["black", "black"], [".", "x"], ["-", "--"], r"$E_{pair}$ BCS vs Exact (N=4)", markersize=0, name="Plots/Epair_N4")
plotEpair([bcs_data_N12, exact_data_N12], ["BCS", "Exact"], ["black", "black"], [".", "x"], ["-", "--"], r"$E_{pair}$ BCS vs Exact (N=12)", markersize=0, name="Plots/Epair_N12")

plotDuration([bcs_data_N4, exact_data_N4], ["BCS", "Exact"], ["black", "black"], [".", "x"], ["-", "--"], r"Duration BCS vs Exact (N=4)", markersize=0, name="Plots/Duration_N4")
plotDuration([bcs_data_N12, exact_data_N12], ["BCS", "Exact"], ["black", "black"], [".", "x"], ["-", "--"], r"Duration BCS vs Exact (N=12)", markersize=0, name="Plots/Duration_N12")

plotOccupationHistogram([bcs_occupation_N4], ["BCS"], ["blue"], r"Occupation numbers BCS (N=4)", name="Plots/Occupation_BCS_N4")
plotOccupationHistogram([exact_occupation_N4], ["Exact"], ["orange"], r"Occupation numbers Exact (N=4)", name="Plots/Occupation_Exact_N4")
plotOccupationHistogram([bcs_occupation_N12], ["BCS"], ["blue"], r"Occupation numbers BCS (N=12)", name="Plots/Occupation_BCS_N12")
plotOccupationHistogram([exact_occupation_N12], ["Exact"], ["orange"], r"Occupation numbers Exact (N=12)", name="Plots/Occupation_Exact_N12")
plotDiffOccupationHistogram([bcs_occupation_N4, exact_occupation_N4], r"Difference in occupation numbers BCS vs Exact (N=4)", name="Plots/Diff_Occupation_N4")
plotDiffOccupationHistogram([bcs_occupation_N12, exact_occupation_N12], r"Difference in occupation numbers BCS vs Exact (N=12)", name="Plots/Diff_Occupation_N12")