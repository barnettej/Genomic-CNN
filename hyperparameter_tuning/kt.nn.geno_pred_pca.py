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

batch = 
epoch = 
train_size = 
valid_size = 
test_size = 
trials = 
repeats = 

initializer = tf.keras.initializers.HeUniform()
# can use usecols to call in only a subset of the snps
annotations = pd.read_csv("/home/user/user/ukbiobank/t2d/data/t2d.annocut.p0.01.train.input.cid.annotations.csv")
annotations_tensor = tf.convert_to_tensor(annotations)
annotations_tensor = tf.cast(annotations_tensor, tf.float32)


train_path = "/home/user/user/ukbiobank/t2d/data/t2d.annocut.p0.01.train.input.with_pca.scale.csv"
valid_path = "/home/user/user/ukbiobank/t2d/data/t2d.annocut.p0.01.valid.input.with_pca.scale.csv"


train_ds = pd.read_csv(train_path)
print(train_ds.head())
x_train, pheno_train, pca_train = train_ds.values[:55168,11:-1], train_ds.values[:55168,10], train_ds.values[:55168,:10]

x_train = tf.cast(tf.convert_to_tensor(x_train), tf.float32)
pca_train = tf.cast(tf.convert_to_tensor(pca_train), tf.float32)
pheno_train = LabelEncoder().fit_transform(pheno_train)

valid_ds = pd.read_csv(valid_path)
x_valid, pheno_valid, pca_valid = valid_ds.values[:11712,11:-1], valid_ds.values[:11712,10], valid_ds.values[:11712,:10]
x_valid = tf.cast(tf.convert_to_tensor(x_valid), tf.float32)
pca_valid = tf.cast(tf.convert_to_tensor(pca_valid), tf.float32)
pheno_valid = LabelEncoder().fit_transform(pheno_valid)



print("x shape:",x_train.shape)
print("pheno shape:", pheno_train.shape)
print("pca shape:", pca_train.shape)

def build_model(hp):
    dense_layers = hp.Int("dense_layers", 1, 3)

    # input genotype data
    input = Input(shape=(x_train.shape[1]), batch_size=batch)
    x = input


    for k in range(dense_layers):
        with hp.conditional_scope("dense_layers", list(range(k + 1, dense_layers + 1))):
            units = hp.Int("dense_" + str(k), 5, 100)
            rate = hp.Float("dense_" + str(k) + "_dropout", min_value=0.0, max_value=0.7)
            x = layers.Dense(units=units, activation="relu",
                             # kernel_regularizer=regularizers.L2(L2)
                             )(x)
            x = layers.Dropout(rate)(x)

    pheno_out = Dense(1, activation="sigmoid", name="pheno_out", kernel_initializer=tf.keras.initializers.GlorotUniform())(x)
    pca_out = Dense(10, name="pca_out", activation="linear")(x)

    losses = {
        "pheno_out": "binary_crossentropy",
        "pca_out": "mse"
    }

    lossWeights = {
        "pheno_out": 0,
        "pca_out": 1
    }

    model = Model(inputs=input, outputs=[pheno_out,pca_out], )

    lr = hp.Float("lr", min_value=1e-6, max_value=1e-3, sampling="log")


    model.compile(optimizer=keras.optimizers.Adam(learning_rate=lr),
                  loss=losses, loss_weights=lossWeights, metrics= {"pheno_out": "AUC", "pca_out": tfa.metrics.RSquare()})


    return (model)


tuner = Hyperband(
    build_model,
    objective=keras_tuner.Objective("val_pca_out_r_square", direction="max"),
    # max_trials=trials,
    # executions_per_trial=repeats,
    project_name="t2d.nn.pca.11.10.22",
    hyperband_iterations= 3,
    # max_epochs = epoch,
    # max_model_size = 500000,
    overwrite=True,
)
# early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', mode="min", patience=5)

# tuner.search_space_summary()


tuner.search(x=x_train, y ={"pheno_out": pheno_train, "pca_out": pca_train},
             validation_data=(x_valid, {"pheno_out": pheno_valid, "pca_out": pca_valid}),
             validation_batch_size=batch,
             validation_steps=(valid_size // batch),
             epochs=epoch,
             batch_size=batch,
             steps_per_epoch=(train_size // batch),
             )

tuner.results_summary()

best_model = tuner.get_best_models(num_models=1)[0]
best_model.summary()
