#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import csv

def main():
#    traj = []
#    with file('../build/real_traj.csv','r') as real_traj:
#        reader = csv.reader(real_traj)
#        for line in reader:
#            traj.append(line)
#        traj = [[float(x) for x in line] for line in traj]
#        plt.figure(100)
#        ax = plt.axes()
#        for state in traj:
#            plt.plot(state[0], state[1], 'o', color='b')
#        for i in range(len(traj)):
#            plt.plot(traj[i][0], traj[i][1], 'o', color='b')
#        plt.plot(traj[0:len(traj)][0],traj[0:len(traj)][1],'o',color='b')
#    plt.show()

    traj = np.loadtxt('../build/real_traj.csv', delimiter=',',dtype="double")
    plt.figure(100)
    plt.plot(traj[:,0], traj[:,1], '-', color='b')
    plt.show()
    traj = np.loadtxt('../build/compare_traj.csv', delimiter=',',dtype="double")
#    plt.figure(101)
#    plt.plot(traj[:,0], traj[:,1], 'ob')
#    plt.plot(traj[:,2], traj[:,3], 'or')
#    plt.show()

if __name__ == '__main__':
    main()
