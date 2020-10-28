from keras import models
from keras import layers
from pandas import read_csv
import random
import numpy as np
import matplotlib.pyplot as plt

data = read_csv('RNA.csv', delimiter=',', header=0)
all_targets = np.array(data["Affinity"])
threshhold = np.median(all_targets)
for i in range(len(all_targets)):
    all_targets[i] = 1 if all_targets[i] > threshhold else 0

all_data = np.zeros((len(data), max([len(data.iloc[i][1]) for i in range(len(data))]), 5))

base_to_index = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
for i in range(len(data)):
    string = str(data.iloc[i][1])
    for j in range(len(string) - 1):
        all_data[i, j, base_to_index[string[j]]] = 1

index = [i for i in range(len(data))]
random.seed(1024)
random.shuffle(index)
all_data = all_data[index]
all_targets = all_targets[index]

boundary = int(len(data)*0.9)
train_data = all_data[:boundary]
train_targets = all_targets[:boundary]
test_data = all_data[boundary:]
test_targets = all_targets[boundary:]


model = models.Sequential()
model.add(layers.Conv1D(256, 5, activation='relu'))
model.add(layers.MaxPooling1D(pool_size=2, strides=2))
model.add(layers.Conv1D(128, 5, activation='relu'))
model.add(layers.GlobalMaxPooling1D())
model.add(layers.Dense(1, activation='sigmoid'))
model.compile(optimizer='rmsprop', loss='binary_crossentropy', metrics=['accuracy'], lr=1e-4)
history = model.fit(train_data, train_targets, epochs=6, batch_size=20, validation_split=0.2)

acc = history.history['accuracy']
val_acc = history.history['val_accuracy']
loss = history.history['loss']
val_loss = history.history['val_loss']

print(model.summary())
epochs = range(1, len(acc) + 1)
plt.plot(epochs, acc, 'bo', label='Training acc')
plt.plot(epochs, val_acc, 'b', label='Validation acc')
plt.title('Training and validation accuracy')
plt.legend()
plt.savefig('acc.png')

plt.figure()
plt.plot(epochs, loss, 'bo', label='Training loss')
plt.plot(epochs, val_loss, 'b', label='Validation loss')
plt.title('Training and validation loss')
plt.legend()
plt.savefig('loss.png')

test_loss, test_acc = model.evaluate(test_data, test_targets)
print(test_acc)