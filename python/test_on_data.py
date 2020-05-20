import pathlib as pt
import pandas as pd
import matplotlib.pyplot as plt
from savgol import ConstantPaddingSavgolFilter

cwd = pt.Path(__file__).parent
resource = pt.Path(cwd/"../resources/BinancBTC_4hr.csv")

assert resource.exists()

data = pd.read_csv(resource)

close = data["Close"].values

plt.plot(close, color="red")

filt1 = ConstantPaddingSavgolFilter(9)
filt2 = ConstantPaddingSavgolFilter(21)

ma1 = filt1(close)

ma2 = filt2(close)

plt.plot(ma1)
plt.plot(ma2)

plt.show()


