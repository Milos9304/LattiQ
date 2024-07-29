import matplotlib.pyplot as plt

#100
#xs=[4,5,6,7,8,9]
#ys_train=[0.545961,0.59289,0.652176,0.625425,0.637182,0.56036]
#ys_test=[0.583482,0.874434,0.734779,0.698994,0.653162,0.54719]

ys_test=[0.802477,0.740654,0.737071,0.706307,0.608728,0.536939,0.494456,0.473054,0.463968,0.462274,0.46494,0.470135,0.476718,0.483995,0.491535,0.499076,0.50646]
xs=[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

fig, ax = plt.subplots()
ax.plot(xs, ys_test, label="Testing performance")
ax.axhline(y=0.5, linestyle='-', color='r', label="Grover's baseline 0.5")
#plt.xlabel("dimensions from 3 to value on x axis")
ax.legend()

# We need to draw the canvas, otherwise the labels won't be positioned and 
# won't have values yet.
fig.canvas.draw()

labels = [item.get_text() for item in ax.get_xticklabels()]
labels = list(map(lambda x: "3-"+x, labels))

ax.set_xticklabels(labels)
plt.title("Pre-trained on 3-9 dimensions, extrapolated. Dimensions in testing vary.")
plt.xlabel("Performance on dataset of dimensions 3-x")
plt.ylabel("alpha")
plt.legend()
plt.show()
