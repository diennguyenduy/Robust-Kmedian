import csv
import random
import pandas as pd
import numpy as np

def random_remove_data(input, output, percentage):
  data = pd.read_csv(input).values

  dimensions = len(data[0])

  def gen_random_numbers_in_range(low, high, n):
      return random.sample(range(low, high), n)

  number = int(len(data)*percentage)
  datum_index = gen_random_numbers_in_range(0, len(data)-1, number)

  with open(output, 'w') as file:
    writer = csv.writer(file)
    header = np.array(['sepal.length','sepal.width','petal.length','petal.width'])
    writer.writerow(header)

    for index, item in enumerate(data):
      number_of_remove_index = random.randint(1,2)
      index_to_remove = random.sample(range(0, 4), number_of_remove_index)
      newDatum = np.zeros(dimensions)
      if index in datum_index:
        for dim in range(dimensions):
          if dim in index_to_remove:
            newDatum[dim] = float(np.nan)
          else:
            newDatum[dim] = item[dim]
      else:
        newDatum = item
      writer.writerow(newDatum)

#Run the function to remove data in the input file
input = 'iris_modified.csv'

percent_arr = [0, 0.1, 0.2, 0.4]

for percent in percent_arr:
  output = 'Missing_{}.csv'.format(int(percent*100))
  random_remove_data(input, output, percent)