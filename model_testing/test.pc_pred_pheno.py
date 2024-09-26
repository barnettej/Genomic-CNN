import os
import tensorflow as tf
import math
from tensorflow import keras
import numpy as np
from tensorflow import data
import keras_tuner
import pandas as pd
from tensorflow.keras import Input, optimizers, regularizers
from tensorflow.keras.layers import Dense, concatenate, Lambda, Multiply, Conv2D, Conv1D, Add, MaxPool1D, Dropout, Reshape, Flatten, \
    MaxPool2D, BatchNormalization, GlobalAveragePooling2D, GlobalAveragePooling1D
from tensorflow.keras import Model
from tensorflow.keras import layers
import keras.backend as K
import datetime
from tensorflow.keras.utils import plot_model
import itertools
from keras_tuner import BayesianOptimization, RandomSearch, Hyperband
from sklearn.preprocessing import LabelEncoder

batch = 64
epoch = 4
train_size = 55168
valid_size = 11712
test_size = 11648


initializer = tf.keras.initializers.HeUniform()


train_path = "/home/user/user/ukbiobank/t2d/data/t2d.annocut.p0.01.train.pca.scale.csv"
valid_path = "/home/user/user/ukbiobank/t2d/data/t2d.annocut.p0.01.valid.pca.scale.csv"
test_path = "/home/user/user/ukbiobank/t2d/data/t2d.annocut.p0.01.test.pca.scale.csv"

train_ds = pd.read_csv(train_path)
print(train_ds.head())
x_train, pheno_train = train_ds.values[:55168,:10], train_ds.values[:55168,-1]
x_train = tf.cast(tf.convert_to_tensor(x_train), tf.float32)
pheno_train = LabelEncoder().fit_transform(pheno_train)

valid_ds = pd.read_csv(valid_path)
x_valid, pheno_valid = valid_ds.values[:11712,:10], valid_ds.values[:11712,-1]
x_valid = tf.cast(tf.convert_to_tensor(x_valid), tf.float32)
pheno_valid = LabelEncoder().fit_transform(pheno_valid)

test_ds = pd.read_csv(test_path)
x_test, pheno_test = test_ds.values[:11648,:10], test_ds.values[:11648,-1]
x_test = tf.cast(tf.convert_to_tensor(x_test), tf.float32)
pheno_test = LabelEncoder().fit_transform(pheno_test)

print("x shape:",x_train.shape)
print("pheno shape:", pheno_train.shape)



# input genotype data
input = Input(shape=(x_train.shape[1]), batch_size=batch)
x = input


x = layers.Dense(units=89, activation = "relu")(x)
x = Dropout(0.10481)(x)
x = layers.Dense(units=38, activation = "relu")(x)
x = Dropout(0.54983)(x)
pheno_out = Dense(1, activation="sigmoid", name="pheno_out", kernel_initializer=tf.keras.initializers.GlorotUniform())(x)

losses = {
    "pheno_out": "binary_crossentropy",
}

lossWeights = {
    "pheno_out": 1.0,
}

model = Model(inputs=input, outputs=[pheno_out], )

lr = 0.00055502


model.compile(optimizer=keras.optimizers.Adam(learning_rate=lr),
              loss=losses, loss_weights=lossWeights, metrics= {"pheno_out": "AUC"})




model.fit(x=x_train, y ={"pheno_out": pheno_train},
             validation_data=(x_valid, {"pheno_out": pheno_valid}),
             validation_batch_size=batch,
             validation_steps=(valid_size // batch),
             epochs=epoch,
             batch_size=batch,
             steps_per_epoch=(train_size // batch),
             )

print("Evaluate on test data")
results = model.evaluate(x=x_test, y={"pheno_out": pheno_test}, batch_size=batch, steps=test_size // batch)
print("test loss, test acc:", results)

print("generate predictions")
predictions = model.predict(x=np.array(x_test), steps=test_size // batch, batch_size = batch)
np.savetxt("test.pc_pred_pheno.csv", predictions,  delimiter=",")

