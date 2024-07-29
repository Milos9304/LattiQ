import matplotlib.pyplot as plt

#1000 angles from paper
#-0.0232639 -0.916444 0.735415 0.550755 0.754262 0.282567 0.865909 0.360061 0.875211 0.284794 1.1415 0.0951161
#xs=[4,5,6,7,8,9,]
#ys_train2=[0.640786,0.684092,0.715216,0.688862,0.642565,0.543488]
#ys_test2=[0.615909,0.609621,0.606176,0.618203,0.57489,0.617128]

#1000 angles from paper
#Final angles: -0.0232639 -0.916444 0.735415 0.550755 0.754262 0.282567 0.865909 0.360061 0.875211 0.284794 1.1415 0.0951161
xs=[4,5,6,7,8,9]
ys_train2=[0.640786,0.684092,0.715216,0.688862,0.642565,0.543488]
ys_test2=[0.695307,0.689919,0.739813,0.694685,0.64252,0.545794]

#1000
xs=[4,5,6,7,8,9]
ys_train=[0.682057,0.688999,0.674238,0.677905,0.656762,0.623537]
ys_test=[0.663139,0.699853,0.686876,0.688862,0.65716,0.624293]

#100
#xs=[4,5,6,7,8,9]
#ys_train=[0.545961,0.59289,0.652176,0.625425,0.637182,0.56036]
#ys_test=[0.583482,0.874434,0.734779,0.698994,0.653162,0.54719]

fig, ax = plt.subplots()
ax.plot(xs, ys_train, label="train 200 instances")
ax.plot(xs, ys_test, label="test 800 instances")
ax.plot(xs, ys_train2, label="train_maxcut_angles 200 inst.")
ax.plot(xs, ys_test2, label="test_maxcut_angles 800 inst.")
ax.axhline(y=0.5, linestyle='-', label="Grover's baseline 0.5")
#plt.xlabel("dimensions from 3 to value on x axis")
ax.legend()

# We need to draw the canvas, otherwise the labels won't be positioned and 
# won't have values yet.
fig.canvas.draw()

labels = [item.get_text() for item in ax.get_xticklabels()]
labels = list(map(lambda x: "3-"+x, labels))

ax.set_xticklabels(labels)

plt.title("Trained and tested on small dimension 3-x. Dimensions in training vary.")
plt.show()
