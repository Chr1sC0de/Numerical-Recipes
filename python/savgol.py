from savgol import SavitzkyGolayFilter
import numpy as np

input_array = np.arange(0,100)

filt = SavitzkyGolayFilter(5, 5, 0, 2)

output = filt(input_array)

print('the code is working')
print(input_array)
print(output)