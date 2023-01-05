import random
import pandas as pd
import numpy as np

def random_remove_data(input, percentage):
  data = pd.read_csv(input).values

  dimensions = len(data[0])

  def gen_random_numbers_in_range(low, high, n):
      return random.sample(range(low, high), n)

  # Select randomly 10 rows of dataset
  number = int(len(data)*percentage)
  datum_index = gen_random_numbers_in_range(0, len(data)-1, number)

  with open('Missing_5.csv', 'w') as file:
    writer = csv.writer(file)
    header = np.array(['x','y','z'])
    writer.writerow(header)

    for index, item in enumerate(data):
        newDatum = np.zeros(dimensions)
        if index in datum_index:
          remove_index = random.randint(0,2)
          for dim in range(dimensions):
            if dim == remove_index:
              newDatum[dim] = float(np.nan)
            else:
              newDatum[dim] = item[dim]
        else:
          newDatum = item
        writer.writerow(newDatum)

input = 'Mall_Customers.csv'

percent_arr = [0, 0.05, 0.1, 0.2]

for percent in percent_arr:
  output = 'Missing_{}.csv'.format(int(percent*100))
  random_remove_data(input, output, percent)