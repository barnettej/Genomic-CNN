import os
import tensorflow as tf
import math
from tensorflow import keras
import tensorflow_probability as tfp
import numpy as np
from tensorflow import data
import keras_tuner
import pandas as pd
from tensorflow.keras import Input, optimizers, regularizers
from tensorflow.keras.layers import Dense, concatenate, Lambda, Multiply, Conv2D, Conv1D, Add, MaxPool1D, Dropout, Reshape, Flatten, \
    MaxPool2D, BatchNormalization, GlobalAveragePooling2D, GlobalAveragePooling1D, AveragePooling1D
from tensorflow.keras import Model
from tensorflow.keras import layers
import keras.backend as K
import datetime
from tensorflow.keras.utils import plot_model
import itertools
from keras_tuner import BayesianOptimization, RandomSearch, Hyperband
from sklearn.preprocessing import LabelEncoder
import tensorflow_addons as tfa


print(tf.__version__)
batch = 64
epoch = 2
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

valid_ds = pd.read_csv(valid_path, dtype=np.float32)
x_valid, pheno_valid, pca_valid = valid_ds.values[:11712,11:], valid_ds.values[:11712,10], valid_ds.values[:11712,:10]
x_valid = tf.cast(tf.convert_to_tensor(x_valid), tf.float32)
pca_valid = tf.cast(tf.convert_to_tensor(pca_valid), tf.float32)
pheno_valid = LabelEncoder().fit_transform(pheno_valid)

test_ds = pd.read_csv(test_path, dtype = np.float32)
x_test, pheno_test, pca_test = test_ds.values[:11648,11:], test_ds.values[:11648,10], test_ds.values[:11648,0:10]
x_test = tf.cast(tf.convert_to_tensor(x_test), tf.float32)
pca_test = tf.cast(tf.convert_to_tensor(pca_test), tf.float32)
pheno_test = LabelEncoder().fit_transform(pheno_test)




print("x shape:",x_train.shape)
print("pheno shape:", pheno_train.shape)
print("pca shape:", pca_train.shape)


class ConvBlock(tf.keras.Model):
    def __init__(self, num_filters, kernel_size, pool_width):
        super(ConvBlock, self).__init__()
        self.conv = Conv1D(num_filters, kernel_size, activation="relu", kernel_initializer=initializer, padding='valid')
        self.batchnorm = BatchNormalization()
        self.avgpool = AveragePooling1D(pool_width)

    def call(self, x):
        x = self.conv(x)
        x = self.batchnorm(x)
        x = self.avgpool(x)

        return x


@tf.custom_gradient
def gradient_reverse(x, lamda=1.0):
    y = tf.identity(x)

    def grad(dy):
        return lamda * -dy, None

    return y, grad

class GradientReversalLayer(tf.keras.layers.Layer):
    def __init__(self, gr_weight):
        super().__init__()
        self.gr_weight = gr_weight

    def call(self, x):
        lamda = self.gr_weight
        return gradient_reverse(x, lamda)


# input genotype data
input = Input(shape=(x_train.shape[1]), batch_size=batch)
geno = tf.expand_dims(input, axis=1)


y = np.ones([input.shape[0], 1, 1]) * annotations_tensor
# switch dimensions of the last 2 dimensions, because conv1d layers create a kernel that uses the last dimension and then runs over the second to last dimension
y = tf.transpose(y, perm=[0, 2, 1])

geno = tf.transpose(geno, perm = [0,2,1])
x = np.ones([1, 1, y.shape[2]]) * geno
x = Multiply()([x, y])
x_geno = np.ones([1, 1, y.shape[2]]) * geno

shared_conv = ConvBlock(66, 17, 10)
x = shared_conv(x)
x_geno = shared_conv(x_geno)
print(x.shape)
shared_conv = ConvBlock(16, 3, 1)
x = shared_conv(x)
x_geno = shared_conv(x_geno)
print(x.shape)



x = Flatten()(x)
x_geno = Flatten()(x_geno)
print("Flattened normal:", x.shape)
print("geno:",x_geno.shape)

shared_dense = layers.Dense(units=98, activation="relu",
                            )
shared_drop = layers.Dropout(0.18176)

x = shared_dense(x)
x_geno = shared_dense(x_geno)
x = shared_drop(x)
x_geno = shared_drop(x_geno)

pheno_out = Dense(1, activation="sigmoid", name="pheno_out", kernel_initializer=tf.keras.initializers.GlorotUniform())(x)
print("pheno out:", pheno_out.shape)
geno_stop = Lambda(lambda x: K.stop_gradient(x))(x_geno)
pheno_out_geno = Dense(1, activation="sigmoid", name="pheno_out_geno", kernel_initializer=tf.keras.initializers.GlorotUniform())(geno_stop)
print("pheno out geno:", pheno_out_geno.shape)


pc_out = Dense(10, name="pc_out", activation="linear")(x)
pc_to_pheno = layers.Dense(units=84, activation="relu")(pc_out)
pc_to_pheno = layers.Dropout(0.24062)(pc_to_pheno)
pc_to_pheno = layers.Dense(units=8, activation="relu")(pc_to_pheno)
pc_to_pheno = layers.Dropout(0.12047)(pc_to_pheno)
pc_to_pheno = layers.Dense(units=90, activation="relu")(pc_to_pheno)
pc_to_pheno = layers.Dropout(0.41068)(pc_to_pheno)
pc_pheno_out = Dense(1, activation="sigmoid", name="pc_pheno_out",
                     kernel_initializer=tf.keras.initializers.GlorotUniform())(pc_to_pheno)

losses = {
    "pheno_out": "binary_crossentropy",
    "pheno_out_geno": "binary_crossentropy",
    "pc_out": "mse",
    "pc_pheno_out": "binary_crossentropy"
}

lossWeights = {
    "pheno_out": 1,
    "pheno_out_geno": 1,
    "pc_out": 1,
    "pc_pheno_out": 1
}

model = Model(inputs=input, outputs=[pheno_out,
                                     pheno_out_geno,
                                     pc_out,
                                     pc_pheno_out
                                     ])

lr = 6.842e-5

model.compile(optimizer=keras.optimizers.Adam(learning_rate=lr),
              loss=losses, loss_weights=lossWeights,
              metrics={"pheno_out": "AUC", "pheno_out_geno": "AUC", "pc_out": tfa.metrics.RSquare(),
                       "pc_pheno_out": "AUC"})

model.fit(x=x_train, y ={"pheno_out": pheno_train, "pheno_out_geno": pheno_train, "pc_out": pca_train, "pc_pheno_out": pheno_train },
             validation_data=(x_valid, {"pheno_out": pheno_valid, "pheno_out_geno": pheno_valid, "pc_out": pca_valid, "pc_pheno_out": pheno_valid}),
             validation_batch_size=batch,
             validation_steps=(valid_size // batch),
             epochs=epoch,
             batch_size=batch,
             steps_per_epoch=(train_size // batch),
             )
model.summary()

print("Evaluate on test data")
results = model.evaluate(x=x_test, y={"pheno_out": pheno_test, "pheno_out_geno": pheno_test, "pc_out": pca_test, "pc_pheno_out": pheno_test}, batch_size=batch, steps=test_size // batch)
print("test loss, test acc:", results)

print("generate predictions")
predictions = model.predict(x=np.array(x_test), steps=test_size // batch, batch_size = batch)
np.savetxt("test.cnn_4_task.csv", predictions[0],  delimiter=",")
np.savetxt("test.cnn_4_task.geno.csv", predictions[1],  delimiter=",")
np.savetxt("test.cnn_4_task.pc_to_pheno.csv", predictions[3],  delimiter=",")
