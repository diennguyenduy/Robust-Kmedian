import numpy as np
from math import sqrt
import pandas as pd
import csv

def CompleteData(k, theta, input, output1, output2):
  data = pd.read_csv(input).values

  dimensions = len(data[0])

  def find_value_and_delta(dim, neighbors, theta):
    max_val = 0
    min_val = 0
    for neighbor in neighbors:
      if np.isnan(neighbor[dim]):
        break
      else:
        max_val = neighbor[dim]
        min_val = neighbor[dim]
        break
      break

    for neighbor in neighbors:
      if neighbor[dim] > max_val:
        max_val = neighbor[dim]
      if neighbor[dim] < min_val:
        min_val = neighbor[dim]

    value = "{:.3f}".format((max_val*(1+theta) + min_val*(1-theta))/2)
    delta = "{:.3f}".format((max_val*(1+theta) - min_val*(1-theta))/2)
    return value, delta

  
  # calculate the Euclidean distance between two vectors
  def euclidean_distance(row1, row2):
    distance = 0.0
    counter1 = 0
    counter2 = 0

    for num_1, num_2 in zip(row1, row2):
      if np.isnan(num_1) or np.isnan(num_2):
        counter1 += 1
      else:
        counter2 += 1
        distance += (num_1-num_2)**2

    return sqrt(distance)/counter2

  # Locate the most similar neighbors
  def get_neighbors(train, test_row, num_neighbors):
    distances = list()
    for train_row in train:
      dist = euclidean_distance(test_row, train_row)
      distances.append((train_row, dist))
      distances.sort(key=lambda tup: tup[1])
    neighbors = list()
    for i in range(num_neighbors):
      neighbors.append(distances[i+1][0]) #use i+1 to skip itself
    return neighbors

  with open(output1, 'w') as file:
    writer = csv.writer(file)
    header = np.array(['x','y','z'])
    writer.writerow(header)
    for datum in data:
      newDatum = np.zeros(dimensions)
      for dim in range(dimensions):
        if np.isnan(datum[dim]):
          neighbors = get_neighbors(data, datum, k)
          result = find_value_and_delta(dim, neighbors, theta)
          newDatum[dim] = result[1]
        else:
          newDatum[dim] = 0
      writer.writerow(newDatum)

  with open(output2, 'w') as file:
    writer = csv.writer(file)
    header = np.array(['x','y','z'])
    writer.writerow(header)
    for datum in data:
      newDatum = np.zeros(dimensions)
      for dim in range(dimensions):
        if np.isnan(datum[dim]):
          neighbors = get_neighbors(data, datum, k)
          result = find_value_and_delta(dim, neighbors, theta)
          newDatum[dim] = result[0]
        else:
          newDatum[dim] = datum[dim]
      writer.writerow(newDatum)

k = 6
percent_arr = [0, 0.05, 0.1, 0.2]
theta_arr = np.array([0, 0.05, 0.10, 0.15])

for percent in percent_arr:
  missing_file = 'Missing_{}.csv'.format(int(percent*100))
  for theta in theta_arr:
    output1 = 'Delta_miss={}_theta={}.csv'.format(int(percent*100), theta)
    output2 = 'CompleteData_miss={}_theta={}.csv'.format(int(percent*100), theta)
    CompleteData(k, theta, missing_file, output1, output2)
