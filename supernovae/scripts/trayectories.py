import numpy as np
import math as m
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import glob
import os
import csv

def main():
	path_data = os.path.join('.','data')
	path_animation = os.path.join(path_data,'animation')
	star_files = []
	fragment_files = []


	for file in os.listdir(path_animation):
		if file.startswith('Frame_star'):
			star_files.append(os.path.join(path_animation,file))
		elif file.startswith('Frame_fragment'):
			fragment_files.append(os.path.join(path_animation,file))

	star = np.zeros((len(star_files)*2,3))
	fragment = np.zeros((len(star_files)*2,3))

	counter = 0
	for file in star_files:
		with open(file) as csv_file:
			csv_reader = csv.reader(csv_file, delimiter=',')
			line_count = 0
			for row in csv_reader:
				if line_count == 0:
				    line_count += 1
				else:
				    star[counter]=(row[0],row[1],row[2])
				    counter +=1
	counter = 0
	for file in fragment_files:
		with open(file) as csv_file:
			csv_reader = csv.reader(csv_file, delimiter=',')
			line_count = 0
			for row in csv_reader:
				if line_count == 0:
				    line_count += 1
				else:
				    fragment[counter]=(row[0],row[1],row[2])
				    counter +=1 

	number_frag = line_count -1

	plot_trayectorie(star, fragment, number_frag)

def plot_trayectorie(stars, fragments, number_frag):
	star1 = np.zeros((int(len(stars)/2),3))
	star2 = np.zeros((int(len(stars)/2),3))

	counter = 0
	for coordinate in stars:
		if counter%2==0 and counter <= len(star1)-1:
			star1[counter]=coordinate
			counter +=1
		elif counter%2!=0 and counter <= len(star2)-1:
			star2[counter]=coordinate
			counter +=1
	
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	ax.scatter(star1[:, 0], star1[:, 1], star1[:, 2], alpha = 0.1)
	ax.scatter(star2[:, 0], star2[:, 1], star2[:, 2], alpha = 0.1)
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	ax.set_xlim3d(-100,100)
	ax.set_ylim3d(-100,100)
	ax.set_zlim3d(-100,100)

	plt.show()

if __name__ == '__main__':
	main()



