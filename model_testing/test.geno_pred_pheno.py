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
import tensorflow_addons as tfa

from keras_tuner import BayesianOptimization, RandomSearch, Hyperband
from sklearn.preprocessing import LabelEncoder

batch = 64
epoch = 4
train_size = 55168
valid_size = 11712
test_size = 11648


initializer = tf.keras.initializers.HeUniform()

annotations = pd.read_csv("/home/user/user/ukbiobank/t2d/data/t2d.annocut.p0.01.train.input.cid.annotations.csv")
annotations_tensor = tf.convert_to_tensor(annotations)
annotations_tensor = tf.cast(annotations_tensor, tf.float32)


train_path = "/home/user/user/ukbiobank/t2d/data/t2d.annocut.p0.01.train.input.with_pca.scale.csv"
valid_path = "/home/user/user/ukbiobank/t2d/data/t2d.annocut.p0.01.valid.input.with_pca.scale.csv"
test_path = "/home/user/user/ukbiobank/t2d/data/t2d.annocut.p0.01.test.input.with_pca.scale.csv"

train_ds = pd.read_csv(train_path, dtype=np.float32)
print(train_ds.head())
x_train, pheno_train, pca_train = train_ds.values[:55168,11:], train_ds.values[:55168,10], train_ds.values[:55168,:10]
x_train = tf.cast(tf.convert_to_tensor(x_train), tf.float32)
pca_train = tf.cast(tf.convert_to_tensor(pca_train), tf.float32)
pheno_train = LabelEncoder().fit_transform(pheno_train)

valid_ds = pd.read_csv(valid_path, dtype = np.float32)
x_valid, pheno_valid, pca_valid = valid_ds.values[:11712,11:], valid_ds.values[:11712,10], valid_ds.values[:11712,:10]
x_valid = tf.cast(tf.convert_to_tensor(x_valid), tf.float32)
pca_valid = tf.cast(tf.convert_to_tensor(pca_valid), tf.float32)
pheno_valid = LabelEncoder().fit_transform(pheno_valid)

test_ds = pd.read_csv(test_path, dtype = np.float32)
x_test, pheno_test, pca_test = test_ds.values[:11648,11:], test_ds.values[:11648,10], test_ds.values[:11648,0:10]
x_test = tf.cast(tf.convert_to_tensor(x_test), tf.float32)
pca_test = tf.cast(tf.convert_to_tensor(pca_test), tf.float32)
pheno_test = LabelEncoder().fit_transform(pheno_test)


@tf.custom_gradient
def gradient_reverse(x, lamda=1.0):
    y = tf.identity(x)

    def grad(dy):
        return lamda * -dy, None

    return y, grad

class GradientReversalLayer(tf.keras.layers.Layer):
    def __init__(self):
        super().__init__()

    def call(selfself, x, lamda=0.00093603):
        return gradient_reverse(x, lamda)



print("x shape:",x_train.shape)
print("pheno shape:", pheno_train.shape)
print("pca shape:", pca_train.shape)



# input genotype data
input = Input(shape=(x_train.shape[1]), batch_size=batch)
x = input
x = layers.Dense(units=77, kernel_regularizer=regularizers.L2(6.7027e-7))(x)
x = layers.Dropout(0.23305)(x)


pheno_out = Dense(1, activation="sigmoid", name="pheno_out", kernel_initializer=tf.keras.initializers.GlorotUniform())(x)
pca_out = Dense(10, name="pca_out", activation="linear")(x)
pca_to_pheno = Lambda(lambda x: K.stop_gradient(x))(pca_out)
pca_to_pheno = layers.Dense(units = 86, activation = "relu")(pca_to_pheno)
pca_to_pheno = layers.Dropout(0.191)(pca_to_pheno)
pca_to_pheno = layers.Dense(units = 38, activation = "relu")(pca_to_pheno)
pca_to_pheno = layers.Dropout(0.24159)(pca_to_pheno)
pca_to_pheno = layers.Dense(units = 9, activation = "relu")(pca_to_pheno)
pca_to_pheno = layers.Dropout(0.16788)(pca_to_pheno)
pc_pheno_out = Dense(1, activation="sigmoid", name="pc_pheno_out", kernel_initializer=tf.keras.initializers.GlorotUniform())(pca_to_pheno)

losses = {
    "pheno_out": "binary_crossentropy",
    "pca_out": "mse",
    "pc_pheno_out": "binary_crossentropy"
}

lossWeights = {
    "pheno_out": 1,
    "pca_out": 0,
    "pc_pheno_out": 0
}

model = Model(inputs=input, outputs=[pheno_out,pca_out, pc_pheno_out], )

lr = 6.3575e-05


model.compile(optimizer=keras.optimizers.Adam(learning_rate=lr),
              loss=losses, loss_weights=lossWeights, metrics= {"pheno_out": "AUC", "pca_out": tfa.metrics.RSquare(), "pc_pheno_out": "AUC"})

model.fit(x=x_train, y={"pheno_out": pheno_train, "pca_out": pca_train, "pc_pheno_out": pheno_train},
          validation_data=(x_valid, {"pheno_out": pheno_valid, "pca_out": pca_valid, "pc_pheno_out": pheno_valid}),
          validation_batch_size=batch,
          validation_steps=(valid_size // batch),
          epochs=epoch,
          batch_size=batch,
          steps_per_epoch=(train_size // batch),
          )

print("Evaluate on test data")
results = model.evaluate(x=x_test, y={"pheno_out": pheno_test, "pca_out": pca_test, "pc_pheno_out": pheno_test}, batch_size=batch, steps=test_size // batch)
print("test loss, test acc:", results)

print("generate predictions")
predictions = model.predict(x=np.array(x_test), steps=test_size // batch, batch_size = batch)
np.savetxt("test.geno_pred_pheno.0_weight.csv", predictions[0],  delimiter=",")
np.savetxt("test.geno_pred_pheno.0_weight_pc_to_pheno.csv", predictions[2], delimiter=",")



